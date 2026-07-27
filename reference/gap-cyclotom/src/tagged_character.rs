//! Character-table workloads over the adaptive tagged cyclotomic kernel.

use crate::libgap::CharacterTable;
use crate::tagged::{Context, Cyclotomic};
use crate::Error;

pub struct PreparedCharacterTable {
    rows: Vec<Vec<Cyclotomic>>,
    weighted_conjugates: Vec<Vec<Cyclotomic>>,
    group_order: i64,
}

pub struct PreparedRepresentationRing {
    rank: usize,
    structure_constants: Vec<i64>,
}

impl PreparedCharacterTable {
    pub fn import(context: &mut Context, table: &CharacterTable) -> Result<Self, Error> {
        let mut rows = Vec::with_capacity(table.rows.len());
        for row in &table.rows {
            let mut imported_row = Vec::with_capacity(row.len());
            for value in row {
                let order = value.order().map_err(|error| Error(error.to_string()))?;
                let coefficients = value
                    .coefficients()
                    .map_err(|error| Error(error.to_string()))?;
                let terms: Vec<_> = coefficients
                    .into_iter()
                    .enumerate()
                    .filter(|(_, (numerator, _))| *numerator != 0)
                    .map(|(exponent, coefficient)| (exponent as u32, coefficient))
                    .collect();
                imported_row.push(context.from_terms(order, &terms)?);
            }
            rows.push(imported_row);
        }

        let mut weighted_conjugates = Vec::with_capacity(rows.len());
        for row in &rows {
            let mut weighted_row = Vec::with_capacity(row.len());
            for (value, class_size) in row.iter().zip(&table.class_sizes) {
                let conjugate = context.conjugate(value)?;
                weighted_row.push(context.scale_integer(&conjugate, *class_size)?);
            }
            weighted_conjugates.push(weighted_row);
        }

        Ok(Self {
            rows,
            weighted_conjugates,
            group_order: table.group_order,
        })
    }

    pub fn tensor_multiplicities(
        &self,
        context: &mut Context,
        lhs: usize,
        rhs: usize,
    ) -> Result<Vec<i64>, Error> {
        let workload = self.rows[lhs].len() * self.weighted_conjugates.len();
        if workload > 25 {
            let products: Vec<_> = self.rows[lhs]
                .iter()
                .zip(&self.rows[rhs])
                .map(|(left, right)| context.mul(left, right))
                .collect::<Result<_, _>>()?;
            return context.sum_product_rows_integer_quotients(
                &products,
                &self.weighted_conjugates,
                self.group_order,
            );
        }

        let mut sums: Vec<Option<Cyclotomic>> = (0..self.rows.len()).map(|_| None).collect();
        for class in 0..self.rows[lhs].len() {
            let product = context.mul(&self.rows[lhs][class], &self.rows[rhs][class])?;
            for (sum, irreducible) in sums.iter_mut().zip(&self.weighted_conjugates) {
                match sum {
                    Some(sum) => {
                        context.mul_add_assign(sum, &product, &irreducible[class])?;
                    }
                    None => {
                        *sum = Some(context.mul(&product, &irreducible[class])?);
                    }
                }
            }
        }
        sums.into_iter()
            .map(|sum| {
                let sum = sum.ok_or_else(|| Error("character table has no classes".into()))?;
                if matches!(&sum, Cyclotomic::Small(_)) {
                    let numerator = integer_value(&sum)?;
                    if numerator % self.group_order == 0 {
                        return Ok(numerator / self.group_order);
                    }
                }
                let multiplicity = context.scale_fraction(&sum, 1, self.group_order)?;
                integer_value(&multiplicity)
            })
            .collect()
    }
}

impl PreparedRepresentationRing {
    pub fn build(table: &PreparedCharacterTable, context: &mut Context) -> Result<Self, Error> {
        let rank = table.rows.len();
        let length = rank
            .checked_mul(rank)
            .and_then(|value| value.checked_mul(rank))
            .ok_or_else(|| Error("representation ring is too large".into()))?;
        let mut structure_constants = vec![0_i64; length];
        for lhs in 0..rank {
            for rhs in lhs..rank {
                let product = table.tensor_multiplicities(context, lhs, rhs)?;
                let forward = (lhs * rank + rhs) * rank;
                structure_constants[forward..forward + rank].copy_from_slice(&product);
                if lhs != rhs {
                    let reverse = (rhs * rank + lhs) * rank;
                    structure_constants[reverse..reverse + rank].copy_from_slice(&product);
                }
            }
        }
        Ok(Self {
            rank,
            structure_constants,
        })
    }

    pub fn rank(&self) -> usize {
        self.rank
    }

    pub fn basis_product(&self, lhs: usize, rhs: usize) -> Result<&[i64], Error> {
        if lhs >= self.rank || rhs >= self.rank {
            return Err(Error(
                "representation-ring basis index is out of range".into(),
            ));
        }
        let start = (lhs * self.rank + rhs) * self.rank;
        Ok(&self.structure_constants[start..start + self.rank])
    }

    pub fn multiply(&self, lhs: &[i64], rhs: &[i64]) -> Result<Vec<i64>, Error> {
        if lhs.len() != self.rank || rhs.len() != self.rank {
            return Err(Error(
                "representation-ring vector has the wrong dimension".into(),
            ));
        }
        let mut result = vec![0_i128; self.rank];
        for (left_index, &left) in lhs.iter().enumerate() {
            if left == 0 {
                continue;
            }
            for (right_index, &right) in rhs.iter().enumerate() {
                if right == 0 {
                    continue;
                }
                let scalar = i128::from(left)
                    .checked_mul(i128::from(right))
                    .ok_or_else(|| Error("representation-ring coefficient overflow".into()))?;
                for (slot, &coefficient) in result
                    .iter_mut()
                    .zip(self.basis_product(left_index, right_index)?)
                {
                    let contribution = scalar
                        .checked_mul(i128::from(coefficient))
                        .ok_or_else(|| Error("representation-ring coefficient overflow".into()))?;
                    *slot = slot
                        .checked_add(contribution)
                        .ok_or_else(|| Error("representation-ring coefficient overflow".into()))?;
                }
            }
        }
        result
            .into_iter()
            .map(|coefficient| {
                i64::try_from(coefficient)
                    .map_err(|_| Error("representation-ring coefficient does not fit i64".into()))
            })
            .collect()
    }
}

fn integer_value(value: &Cyclotomic) -> Result<i64, Error> {
    if value.order() != 1 {
        return Err(Error("character multiplicity is not rational".into()));
    }
    let terms = value.rational_terms();
    match terms.as_slice() {
        [] => Ok(0),
        [(0, coefficient)] if coefficient.denom() == &1 => coefficient
            .numer()
            .to_i64()
            .ok_or_else(|| Error("character multiplicity does not fit i64".into())),
        _ => Err(Error("character multiplicity is not an integer".into())),
    }
}
