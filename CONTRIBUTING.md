# Contributing

FASTA Lens is a new project. Contributions and independent reviews are welcome.

1. Open an issue describing the problem or proposed behavior. For bugs, include
   your Python version, operating system, command, and a small synthetic input.
2. Fork the repository and make a branch for one focused change.
3. Add a regression test for a bug or behavior change. Keep the runtime dependency free.
4. Run `python -m unittest discover -s tests -v` and update docs where behavior changes.
5. Submit a pull request explaining the behavior, rationale, and validation.

Do not upload private, patient, proprietary, or unpublished sequence data.
Use synthetic sequences in tests and public examples.

Report JSON changes should preserve existing fields or explicitly change the
schema version. Keep metric denominators and threshold boundaries documented.

AI-assisted contributions are welcome. Explain substantial generated changes,
review them, and include meaningful verification. Do not invent test outcomes,
benchmarks, users, or biological validation.

Be respectful and focus reviews on code and evidence.
