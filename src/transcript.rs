use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript;

pub fn append<Item>(transcript: &mut Transcript, label: &'static [u8], item: Item)
where
    Item: CanonicalSerialize
{
    let mut bytes = Vec::new();
    item.serialize(&mut bytes).unwrap();
    transcript.append_message(label, &bytes)
}

pub fn challenge_scalar<F>(transcript: &mut Transcript, label: &'static [u8]) -> F
where
    F: PrimeField
{
    let size = F::size_in_bits() / 8;
    let mut buf = vec![0u8; size];
    transcript.challenge_bytes(label, &mut buf);
    F::from_random_bytes(&buf).unwrap()
}