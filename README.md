# caulk
Implementation of the Caulk protocol: https://eprint.iacr.org/2022/621.pdf

## Status
Currently only the "blind evaluation protocol" for a single opening is implemented.

## To do
    - verify pedersen commitment opening as in Section 4.7
    - verify that z_comm is correctly constructed as in Section 6.2

## Running
To run the blind evalutation protocol for a single opening, run:
`cargo test single`
