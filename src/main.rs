#![allow(non_snake_case)]
use std::cmp::max;

use ark_bls12_381::{Bls12_381, Fr, G1Affine, G1Projective, G2Affine};
use ark_ec::{bls12::Bls12, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::*;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_poly_commit::kzg10::{Powers, KZG10};

mod tools;
use crate::tools::*;

type UniPoly381 = DensePolynomial<<Bls12_381 as PairingEngine>::Fr>;
type KzgBls12_381 = KZG10<Bls12_381, UniPoly381>;

pub struct Proof {
    // Prover's Pedersen commitment to their value
    cm: G1Affine,

    z_comm: G2Affine,
    T_comm: G1Affine,
    S_comm: G2Affine,
}

pub struct Setup {
    // multiplicative subgroup H
    domain: GeneralEvaluationDomain<Fr>,

    // generator of Group 1
    g1: G1Affine,

    // generator of Group 2
    g2: G2Affine,

    // generator of Group 2 times the secret used in KZG setup
    xg2: G2Affine,

    // common point used in Pedersen Hash. For randomness `r`, Ped_r(v) = v*[g1] + r*[h]
    h: G1Affine,

    // KZG SRS: {g1, x*g1, x^2*g1, ..., x^{2^k}*g1}
    ck: Powers<'static, Bls12<ark_bls12_381::Parameters>>,
}

pub fn setup(n: usize) -> Setup {
    let rng = &mut rand::thread_rng();

    let domain: GeneralEvaluationDomain<Fr> = GeneralEvaluationDomain::new(n + 2).unwrap();

    let ck_size = max(n, 2 * domain.size() + 3);

    let srs = KzgBls12_381::setup(max(n, ck_size), true, rng).unwrap();

    // trim down to size.
    let (ck, vk) = trim::<Bls12_381, UniPoly381>(srs, ck_size);

    let g1: G1Affine = ck.powers_of_g[0];
    let g2: G2Affine = vk.h;
    let xg2: G2Affine = vk.beta_h;
    let h: G1Affine = G1Projective::rand(rng).into_affine();

    Setup {
        domain,
        g1,
        g2,
        xg2,
        h,
        ck,
    }
}

fn prove(setup: &Setup, entries: &[Fr], index: usize, value: Fr) -> Proof {
    let rng = &mut rand::thread_rng();

    let domain = setup.domain;
    let g1 = setup.g1;
    let g2 = setup.g2;
    let xg2 = setup.xg2;
    let h = setup.h;
    let ck = &setup.ck;

    // Prover's computes `cm` which is a Pedersen commitment to `value` using randomness `r`
    //  cm = value*[g1] + r*[h]
    let r = Fr::rand(rng);
    let cm = (g1.mul(value) + h.mul(r)).into_affine();

    // prover chooses random `a` and `s`
    let a = Fr::rand(rng);
    let s = Fr::rand(rng);

    // Lagrange basis polynomial encoding the entries
    let C_poly = DensePolynomial::from_coefficients_vec(domain.ifft(entries));

    // prover computes q, T, and S
    let v_poly = DensePolynomial::from_coefficients_vec(vec![value]);
    let d_poly = DensePolynomial::from_coefficients_vec(vec![-domain.element(index), Fr::one()]);
    let q_poly = &(&C_poly - &v_poly) / &d_poly;

    // Prover computes commitment to `q`
    let q_comm = commit(&q_poly, ck);

    // Prover computes commitment to `T` in group 1
    // [T]_1 = [a^(-1)*q + s*h]_1 = a^(-1)*[q]_1 + s*[h]_1
    let T_comm = (q_comm.mul(a.inverse().unwrap()) + h.mul(s)).into_affine();

    // Prover computes commitment to `z` in group 2
    // [z]_2 = [a(X-omega^i)]_2 = a*[X]_2 - a*omega^i[1]_2
    let z_comm = (xg2.mul(a) - g2.mul(a * domain.element(index))).into_affine();

    // Prover computes commitment to `S` in group 2
    // [S]_2 = [-r - s*z]_2
    let S_comm = (g2.mul(-r) - z_comm.mul(s)).into_affine();

    Proof {
        cm,
        z_comm,
        T_comm,
        S_comm,
    }
}

pub fn verify(setup: &Setup, entries: &[Fr], proof: Proof) -> bool {
    // Lagrange basis polynomial encoding the entries
    let C_poly = DensePolynomial::from_coefficients_vec(setup.domain.ifft(entries));
    let C_comm = commit(&C_poly, &setup.ck).into_projective();

    let cm = proof.cm.into_projective();
    let T_comm = proof.T_comm;
    let z_comm = proof.z_comm;
    let S_comm = proof.S_comm;

    let g2 = setup.g2;
    let h = setup.h;

    let pairing_left = Bls12_381::pairing(C_comm - cm, g2);
    let pairing_right = Bls12_381::pairing(T_comm, z_comm) * Bls12_381::pairing(h, S_comm);

    pairing_left == pairing_right
}

fn main() {
    let rng = &mut rand::thread_rng();

    let num_of_entries = 100;

    // we have a list of __field_elements__ representing notes
    let entries: Vec<Fr> = (0..num_of_entries)
        .into_iter()
        .map(|_i| Fr::rand(rng))
        .collect();

    let setup = setup(num_of_entries);

    // the Prover owns the entry at index 47
    let prover_entry_index = 47usize;
    let prover_entry = entries[prover_entry_index];

    // Prover constructs a proof which contains a hiding Pedersen commitment to their entry
    let proof = prove(&setup, &entries, prover_entry_index, prover_entry);

    // Verifier verifies the Prover's Pedersen commitment corresponds to some entry in the list
    let result = verify(&setup, &entries, proof);

    //  NOTE: as currently written the verifier code is incomplete
    //  TODO:
    //      - verify pedersen commitment opening as in Section 4.7
    //      - verify that z_comm is correctly constructed as in Section 6.2

    if result {
        println!("Proof verified!")
    } else {
        println!("Proof invalid")
    }
}
