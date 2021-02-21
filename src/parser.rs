use crate::fields::GenericCyclotomic;
use crate::fields::Z;
use nom::sequence::delimited;
use nom::IResult;
use nom::{
    branch::alt,
    character::complete::{char, digit1},
    combinator::{map, map_res, opt},
    sequence::tuple,
};
use std::collections::HashMap;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
enum Expr {
    Root(u64),
    Integer(i64),
    Rational(i64, u64),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mult(Box<Expr>, Box<Expr>),
    Power(Box<Expr>, u64),
    List(Vec<Expr>),
}

// TODO: this is kind of hard to read...
fn integer(s: &str) -> IResult<&str, i64> {
    map_res(tuple((opt(char('-')), digit1)), |(sign, x)| {
        FromStr::from_str((sign.map_or("", |s| "-").to_owned() + x).as_str())
    })(s)
}

#[test]
fn test_integer() {
    assert_eq!(integer("123"), Ok(("", 123)));
    assert_eq!(integer("-1"), Ok(("", -1)));
}

fn fraction(s1: &str) -> IResult<&str, Expr> {
    let (s2, (numer, _, denom)) = tuple((integer, char('/'), integer))(s1)?;
    Ok((s2, Expr::Rational(numer, denom as u64)))
}

#[test]
fn test_fraction() {
    assert_eq!(fraction("123/456"), Ok(("", Expr::Rational(123, 456))));
    assert_eq!(fraction("-123/456"), Ok(("", Expr::Rational(-123, 456))));
    assert_eq!(fraction("1/0"), Ok(("", Expr::Rational(1, 0))));
}

fn coeff(s1: &str) -> IResult<&str, Expr> {
    let (s2, q) = opt(fraction)(s1)?;

    match q {
        None => map(integer, Expr::Integer)(s1),
        Some(rat) => Ok((s2, rat)),
    }
}

#[test]
fn test_coeff() {
    assert_eq!(coeff("123"), Ok(("", Expr::Integer(123))));
    assert_eq!(coeff("-1"), Ok(("", Expr::Integer(-1))));
    assert_eq!(coeff("123/456"), Ok(("", Expr::Rational(123, 456))));
    assert_eq!(coeff("-123/456"), Ok(("", Expr::Rational(-123, 456))));
    assert_eq!(coeff("1/0"), Ok(("", Expr::Rational(1, 0))));
}

fn root_no_power(s1: &str) -> IResult<&str, Expr> {
    let (s2, (_, k)) = tuple((char('E'), delimited(char('('), integer, char(')'))))(s1)?;
    Ok((s2, Expr::Root(k as u64)))
}

#[test]
fn test_root_no_power() {
    assert_eq!(root_no_power("E(6)"), Ok(("", Expr::Root(6))));
    assert_eq!(root_no_power("E(0)"), Ok(("", Expr::Root(0))));
    assert_eq!(root_no_power("E(1234)"), Ok(("", Expr::Root(1234))));
}

fn root_power(s1: &str) -> IResult<&str, Expr> {
    let (s2, base) = root_no_power(s1)?;
    let (s3, (_, power)) = tuple((char('^'), integer))(s2)?;
    Ok((s3, Expr::Power(Box::new(base), power as u64)))
}

#[test]
fn test_root_power() {
    assert_eq!(
        root_power("E(6)^2"),
        Ok(("", Expr::Power(Box::new(Expr::Root(6)), 2)))
    );
    assert_eq!(
        root_power("E(123)^56"),
        Ok(("", Expr::Power(Box::new(Expr::Root(123)), 56)))
    );
    assert_eq!(
        root_power("E(0)^0"),
        Ok(("", Expr::Power(Box::new(Expr::Root(0)), 0)))
    );
}

fn root(s1: &str) -> IResult<&str, Expr> {
    alt((root_power, root_no_power))(s1)
}

#[test]
fn test_root() {
    assert_eq!(
        root("E(6)^2"),
        Ok(("", Expr::Power(Box::new(Expr::Root(6)), 2)))
    );
    assert_eq!(
        root("E(123)^56"),
        Ok(("", Expr::Power(Box::new(Expr::Root(123)), 56)))
    );
    assert_eq!(
        root("E(0)^0"),
        Ok(("", Expr::Power(Box::new(Expr::Root(0)), 0)))
    );
    assert_eq!(root("E(6)"), Ok(("", Expr::Root(6))));
    assert_eq!(root("E(0)"), Ok(("", Expr::Root(0))));
    assert_eq!(root("E(1234)"), Ok(("", Expr::Root(1234))));
}

fn rat_root(s1: &str) -> IResult<&str, Expr> {
    let (s2, (coeff, _, root)) = tuple((coeff, char('*'), root))(s1)?;
    Ok((s2, Expr::Mult(Box::new(coeff), Box::new(root))))
}

fn coeff_root(s1: &str) -> IResult<&str, Expr> {
    let (s2, maybe_rat_root) = opt(rat_root)(s1)?;
    match maybe_rat_root {
        Some(rat_root) => Ok((s2, rat_root)),
        None => {
            let (s3, (_, root)) = tuple((char('-'), root))(s1)?;
            Ok((s3, Expr::Mult(Box::new(Expr::Integer(-1)), Box::new(root))))
        }
    }
}

#[test]
fn test_coeff_root() {
    assert_eq!(
        coeff_root("12*E(7)"),
        Ok((
            "",
            Expr::Mult(Box::new(Expr::Integer(12)), Box::new(Expr::Root(7)))
        ))
    );
    assert_eq!(
        coeff_root("1*E(6)^2"),
        Ok((
            "",
            Expr::Mult(
                Box::new(Expr::Integer(1)),
                Box::new(Expr::Power(Box::new(Expr::Root(6)), 2))
            )
        ))
    );
    assert_eq!(
        coeff_root("2/40*E(16)^12"),
        Ok((
            "",
            Expr::Mult(
                Box::new(Expr::Rational(2, 40)),
                Box::new(Expr::Power(Box::new(Expr::Root(16)), 12))
            )
        ))
    );
    assert_eq!(
        coeff_root("-E(60)^29"),
        Ok((
            "",
            Expr::Mult(
                Box::new(Expr::Integer(-1)),
                Box::new(Expr::Power(Box::new(Expr::Root(60)), 29))
            )
        ))
    );
}

pub fn parse_element(s: &str) -> Option<GenericCyclotomic> {
    Some(GenericCyclotomic {
        exp_coeffs: vec![(Z::from(0), (123, 1))].into_iter().collect(),
        order: Z::from(1),
    })
}

pub fn parse_vector(s: &str) -> Option<Vec<GenericCyclotomic>> {
    Some(vec![GenericCyclotomic {
        exp_coeffs: HashMap::new(),
        order: Z::from(1),
    }])
}

pub fn parse_matrix(s: &str) -> Option<Vec<Vec<GenericCyclotomic>>> {
    Some(vec![vec![GenericCyclotomic {
        exp_coeffs: HashMap::new(),
        order: Z::from(1),
    }]])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_integer() {
        assert_eq!(
            parse_element("123").unwrap(),
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (123, 1))].into_iter().collect(),
                order: Z::from(1)
            }
        )
    }

    #[test]
    fn parses_rational() {
        assert_eq!(
            parse_element("12/4").unwrap(),
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (12, 4))].into_iter().collect(),
                order: Z::from(1)
            }
        )
    }
}
