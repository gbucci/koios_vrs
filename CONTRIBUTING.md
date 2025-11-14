# Contributing to KOIOS-VRS

Thank you for your interest in contributing to KOIOS-VRS (Kinship Open Inference Overhaul System - Variant Reference Selection)!

## Branching Strategy

We use a protected branching strategy to maintain stable releases on Zenodo while enabling active development:

### Branch Structure

- **`main`** - Protected branch containing stable, versioned releases archived on Zenodo
  - Only maintainers can merge to this branch
  - Each commit represents a versioned release
  - DO NOT push directly to this branch

- **`develop`** - Main development branch
  - All feature branches should be created from and merged into `develop`
  - Collaborative work happens here
  - Periodically merged into `main` for releases

- **`feature/*`** - Feature branches for specific developments
  - Created from `develop`
  - Merged back into `develop` via pull request
  - Naming convention: `feature/descriptive-name`

## How to Contribute

### 1. Clone the Repository

```bash
git clone https://github.com/gbucci/koios_vrs.git
cd koios_vrs
```

### 2. Create a Feature Branch

Always create your feature branch from `develop`:

```bash
git checkout develop
git pull origin develop
git checkout -b feature/your-feature-name
```

### 3. Make Your Changes

- Write clear, commented code
- Follow the existing code style
- Test your changes thoroughly
- Update documentation as needed

### 4. Commit Your Changes

Write clear, descriptive commit messages:

```bash
git add .
git commit -m "Brief description of changes"
```

### 5. Push Your Branch

```bash
git push origin feature/your-feature-name
```

### 6. Create a Pull Request

- Go to the GitHub repository
- Create a pull request from your feature branch to `develop`
- Provide a clear description of your changes
- Reference any related issues

## Code Guidelines

### R Code Style

- Use meaningful variable and function names
- Comment complex logic
- Follow tidyverse style guide where applicable
- Include roxygen2 documentation for functions

### Testing

- Test your code with different VCF inputs
- Verify that existing functionality still works
- Include example usage in your PR description

### Documentation

- Update README.md if adding new features
- Add comments to explain complex algorithms
- Update function documentation

## Reporting Issues

If you find a bug or have a suggestion:

1. Check if the issue already exists
2. Create a new issue with:
   - Clear title
   - Detailed description
   - Steps to reproduce (for bugs)
   - Expected vs actual behavior
   - System information (OS, R version)

## Questions?

Feel free to open an issue for questions about contributing or using KOIOS-VRS.

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).
