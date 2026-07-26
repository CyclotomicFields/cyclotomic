//! Character-table workloads over the adaptive tagged cyclotomic kernel.

use crate::libgap::CharacterTable;
use crate::tagged::{Context, Cyclotomic};
use crate::Error;

pub struct PreparedCharacterTable {
    rows: Vec<Vec<Cyclotomic>>,
    weighted_conjugates: Vec<Vec<Cyclotomic>>,
    group_order: i64,
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
        let mut sums: Vec<Option<Cyclotomic>> = (0..self.rows.len()).map(|_| None).collect();

        for class in 0..self.rows[lhs].len() {
            let product = context.mul(&self.rows[lhs][class], &self.rows[rhs][class])?;
            for (sum, irreducible) in sums.iter_mut().zip(&self.weighted_conjugates) {
                let term = context.mul(&product, &irreducible[class])?;
                *sum = Some(match sum.take() {
                    Some(mut sum) => {
                        context.add_assign(&mut sum, &term)?;
                        sum
                    }
                    None => term,
                });
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
