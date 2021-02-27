use crate::fields::GenericCyclotomic;
use crate::fields::Z;
use nom::multi::{fold_many0, separated_list0};
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

#[derive(Debug, PartialEq, Eq, Clone)]
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

fn term(s1: &str) -> IResult<&str, Expr> {
    alt((root, coeff_root, coeff))(s1)
}

#[test]
fn test_term() {
    assert_eq!(
        term("-E(60)^29"),
        Ok((
            "",
            Expr::Mult(
                Box::new(Expr::Integer(-1)),
                Box::new(Expr::Power(Box::new(Expr::Root(60)), 29))
            )
        ))
    );
    assert_eq!(
        term("E(60)^29"),
        Ok(("", Expr::Power(Box::new(Expr::Root(60)), 29)))
    );
}

fn expr(s1: &str) -> IResult<&str, Expr> {
    let (s2, leading_term) = term(s1)?;

    fold_many0(
        tuple((alt((char('+'), char('-'))), term)),
        leading_term,
        |acc, (op, next_term)| {
            if op == '+' {
                Expr::Add(Box::new(acc), Box::new(next_term))
            } else {
                Expr::Sub(Box::new(acc), Box::new(next_term))
            }
        },
    )(s2)
}

#[test]
fn test_expr() {
    assert_eq!(
        expr("-E(60)^29"),
        Ok((
            "",
            Expr::Mult(
                Box::new(Expr::Integer(-1)),
                Box::new(Expr::Power(Box::new(Expr::Root(60)), 29))
            )
        ))
    );
    assert_eq!(
        expr("E(6)^2+E(6)^3"),
        Ok((
            "",
            Expr::Add(
                Box::new(Expr::Power(Box::new(Expr::Root(6)), 2)),
                Box::new(Expr::Power(Box::new(Expr::Root(6)), 3))
            )
        ))
    );
    assert_eq!(
        expr("-E(5)^2+1/2*E(5)^3-E(5)"),
        Ok((
            "",
            Expr::Sub(
                Box::new(Expr::Add(
                    Box::new(Expr::Mult(
                        Box::new(Expr::Integer(-1)),
                        Box::new(Expr::Power(Box::new(Expr::Root(5)), 2))
                    )),
                    Box::new(Expr::Mult(
                        Box::new(Expr::Rational(1, 2)),
                        Box::new(Expr::Power(Box::new(Expr::Root(5)), 3))
                    ))
                )),
                Box::new(Expr::Root(5))
            )
        ))
    );
}

fn extract_order_power(root: Box<Expr>) -> (u64, u64) {
    match *root {
        Expr::Root(k) => (k, 1),
        Expr::Power(root, power) => (extract_order_power(root).0, power),
        _ => panic!("not a root!"),
    }
}

fn consume_term(term: Box<Expr>, sign: i64, result: &mut GenericCyclotomic) {
    match *term {
        Expr::Root(_) | Expr::Power(_, _) => {
            let (order, power) = extract_order_power(term);
            result.order = Z::from(order);
            result.exp_coeffs.insert(Z::from(power), (sign, 1));
        }
        Expr::Integer(x) => {
            result.exp_coeffs.insert(Z::from(0), (sign * x, 1));
        }
        Expr::Rational(p, q) => {
            result.exp_coeffs.insert(Z::from(0), (sign * p, q));
        }
        Expr::Mult(coeff, root) => {
            let rational = match *coeff {
                Expr::Integer(x) => (sign * x, 1),
                Expr::Rational(p, q) => (sign * p, q),
                _ => (1, 1),
            };
            let (order, power) = extract_order_power(root);
            result.order = Z::from(order);
            result.exp_coeffs.insert(Z::from(power), rational);
        }
        _ => {}
    }
}

fn consume_terms(expression: Box<Expr>, result: &mut GenericCyclotomic) {
    match *expression {
        Expr::Add(left, right) => {
            consume_term(right, 1, result);
            consume_terms(left, result)
        }
        Expr::Sub(left, right) => {
            consume_term(right, -1, result);
            consume_terms(left, result)
        }
        _ => consume_term(expression, 1, result),
    }
}

fn expr2cyc(expression: Expr) -> Option<GenericCyclotomic> {
    let mut result = GenericCyclotomic {
        exp_coeffs: HashMap::new(),
        order: Z::from(1),
    };

    consume_terms(Box::new(expression), &mut result);
    Some(result)
}

pub fn parse_element(s: &str) -> Option<GenericCyclotomic> {
    let parsed = expr(s);

    match parsed {
        Ok((_, expression)) => expr2cyc(expression),
        Err(_) => None,
    }
}

#[test]
fn parses_single_roots() {
    assert_eq!(
        parse_element("123").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(0), (123, 1))].into_iter().collect(),
            order: Z::from(1)
        }
    );
    assert_eq!(
        parse_element("12/4").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(0), (12, 4))].into_iter().collect(),
            order: Z::from(1)
        }
    );
    assert_eq!(
        parse_element("E(5)").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(1), (1, 1))].into_iter().collect(),
            order: Z::from(5)
        }
    );
    assert_eq!(
        parse_element("E(6)^2").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(2), (1, 1))].into_iter().collect(),
            order: Z::from(6)
        }
    );
    assert_eq!(
        parse_element("-E(5)").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(1), (-1, 1))].into_iter().collect(),
            order: Z::from(5)
        }
    );
    assert_eq!(
        parse_element("-E(7)^3").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(3), (-1, 1))].into_iter().collect(),
            order: Z::from(7)
        }
    );
    assert_eq!(
        parse_element("2*E(7)^3").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(3), (2, 1))].into_iter().collect(),
            order: Z::from(7)
        }
    );
    assert_eq!(
        parse_element("2/3*E(700)^345").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(345), (2, 3))].into_iter().collect(),
            order: Z::from(700)
        }
    );
}

#[test]
fn parses_multiple_roots() {
    assert_eq!(
        parse_element("E(3)+E(3)^2").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(1), (1, 1)), (Z::from(2), (1, 1))]
                .into_iter()
                .collect(),
            order: Z::from(3)
        }
    );
    assert_eq!(
        parse_element("-E(3)-E(3)^2").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(1), (-1, 1)), (Z::from(2), (-1, 1))]
                .into_iter()
                .collect(),
            order: Z::from(3)
        }
    );
    assert_eq!(
        parse_element("1/2*E(4)-2/3*E(4)^2").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![(Z::from(1), (1, 2)), (Z::from(2), (-2, 3))]
                .into_iter()
                .collect(),
            order: Z::from(4)
        }
    );
    assert_eq!(
        parse_element("-1/2*E(5)-2/30*E(5)^2+4/5*E(5)^3").unwrap(),
        GenericCyclotomic {
            exp_coeffs: vec![
                (Z::from(1), (-1, 2)),
                (Z::from(2), (-2, 30)),
                (Z::from(3), (4, 5))
            ]
            .into_iter()
            .collect(),
            order: Z::from(5)
        }
    );
}

fn vector(s: &str) -> IResult<&str, Vec<Expr>> {
    delimited(char('['), separated_list0(char(','), expr), char(']'))(s)
}

pub fn parse_vector(s: &str) -> Option<Vec<GenericCyclotomic>> {
    match vector(s) {
        Ok((_, cycs)) => cycs.into_iter().map(expr2cyc).collect(),
        Err(_) => None,
    }
}

#[test]
fn test_vector() {
    assert_eq!(
        parse_vector("[1,2,3]").unwrap(),
        vec![
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (1, 1))].into_iter().collect(),
                order: Z::from(1)
            },
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (2, 1))].into_iter().collect(),
                order: Z::from(1)
            },
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (3, 1))].into_iter().collect(),
                order: Z::from(1)
            },
        ]
    );
    assert_eq!(
        parse_vector("[E(3)+E(3)^2,2,-1/2*E(5)-2/30*E(5)^2+4/5*E(5)^3]").unwrap(),
        vec![
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(1), (1, 1)), (Z::from(2), (1, 1))]
                    .into_iter()
                    .collect(),
                order: Z::from(3)
            },
            GenericCyclotomic {
                exp_coeffs: vec![(Z::from(0), (2, 1))].into_iter().collect(),
                order: Z::from(1)
            },
            GenericCyclotomic {
                exp_coeffs: vec![
                    (Z::from(1), (-1, 2)),
                    (Z::from(2), (-2, 30)),
                    (Z::from(3), (4, 5))
                ]
                .into_iter()
                .collect(),
                order: Z::from(5)
            }
        ]
    );
}

fn matrix(s: &str) -> IResult<&str, Vec<Vec<Expr>>> {
    delimited(char('['), separated_list0(char(','), vector), char(']'))(s)
}

pub fn parse_matrix(s: &str) -> Option<Vec<Vec<GenericCyclotomic>>> {
    match matrix(s) {
        Ok((_, vecs)) => vecs
            .into_iter()
            .map(|v| v.into_iter().map(expr2cyc).collect())
            .collect(),
        Err(_) => None,
    }
}

#[test]
fn test_matrix() {
    let one = GenericCyclotomic {
        exp_coeffs: vec![(Z::from(0), (1, 1))].into_iter().collect(),
        order: Z::from(1),
    };
    let zero = GenericCyclotomic {
        exp_coeffs: vec![(Z::from(0), (0, 1))].into_iter().collect(),
        order: Z::from(1),
    };
    assert_eq!(
        parse_matrix("[[1,0,0],[0,1,0],[0,0,1]]").unwrap(),
        vec![
            vec![one.clone(), zero.clone(), zero.clone()],
            vec![zero.clone(), one.clone(), zero.clone()],
            vec![zero.clone(), zero.clone(), one.clone()]
        ]
    )
}
