#![allow(dead_code)]
use ark_bls12_381::{Bls12_381, Fr, G1Affine};
use ark_ec::{bls12::Bls12, PairingEngine};
use ark_poly::{univariate::DensePolynomial, GeneralEvaluationDomain, UVPolynomial};
use ark_poly_commit::kzg10::*;

type UniPoly381 = DensePolynomial<<Bls12_381 as PairingEngine>::Fr>;
type KzgBls12_381 = KZG10<Bls12_381, UniPoly381>;

// Reduces full srs down to smaller srs for smaller polynomials
// Copied from arkworks library (where same function is private)
pub fn trim<E: PairingEngine, P: UVPolynomial<E::Fr>>(
    srs: UniversalParams<E>,
    mut supported_degree: usize,
) -> (Powers<'static, E>, VerifierKey<E>) {
    if supported_degree == 1 {
        supported_degree += 1;
    }
    let pp = srs;
    let powers_of_g = pp.powers_of_g[..=supported_degree].to_vec();
    let powers_of_gamma_g = (0..=supported_degree)
        .map(|i| pp.powers_of_gamma_g[&i])
        .collect();

    let powers = Powers {
        powers_of_g: ark_std::borrow::Cow::Owned(powers_of_g),
        powers_of_gamma_g: ark_std::borrow::Cow::Owned(powers_of_gamma_g),
    };
    let vk = VerifierKey {
        g: pp.powers_of_g[0],
        gamma_g: pp.powers_of_gamma_g[&0],
        h: pp.h,
        beta_h: pp.beta_h,
        prepared_h: pp.prepared_h.clone(),
        prepared_beta_h: pp.prepared_beta_h.clone(),
    };
    (powers, vk)
}

pub fn commit(
    poly: &DensePolynomial<Fr>,
    ck: &Powers<Bls12<ark_bls12_381::Parameters>>,
) -> G1Affine {
    let (comm, _) = KzgBls12_381::commit(ck, poly, None, None).unwrap();
    comm.0
}

pub fn format_poly_coeffs(poly: &DensePolynomial<Fr>) -> String {
    poly.coeffs
        .iter()
        .enumerate()
        .map(|(i, c)| format!("{} * x^{}", c, i))
        .collect::<Vec<String>>()
        .join("\n")
}

pub fn format_poly_evals(
    poly: &DensePolynomial<Fr>,
    domain: GeneralEvaluationDomain<Fr>,
) -> String {
    poly.clone()
        .evaluate_over_domain(domain)
        .evals
        .iter()
        .map(|c| format!("{}", c))
        .collect::<Vec<String>>()
        .join("\n")
}
