# CI/CD Integration Guide

This document describes how to integrate the code quality checks into your CI/CD pipeline.

## Overview

The code quality check scripts support CI mode, which:
- Runs non-interactively (no user prompts)
- Exits with appropriate error codes
- Can be integrated into any CI/CD system

## Exit Codes

The `run_all_checks.py` script uses the following exit codes:

- `0`: All checks passed
- `1-5`: Number of failed tools (1 = 1 tool failed, 2 = 2 tools failed, etc.)

## GitHub Actions

### Basic Workflow

Create `.github/workflows/code-quality.yml`:

```yaml
name: Code Quality Checks

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]

jobs:
  quality:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      
      - name: Install dependencies
        run: |
          pip install -r requirements.txt
          pip install -r requirements-dev.txt
      
      - name: Run code quality checks
        run: |
          python scripts/run_all_checks.py --ci
      
      - name: Upload reports
        if: always()
        uses: actions/upload-artifact@v3
        with:
          name: lint-reports
          path: lint_reports/
```

### Advanced Workflow with Separate Jobs

```yaml
name: Code Quality Checks

on:
  push:
    branches: [ main, develop ]
  pull_request:
    branches: [ main, develop ]

jobs:
  formatting:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install black
      - run: black --check .

  style:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install flake8
      - run: flake8 compass services

  tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install -r requirements.txt
      - run: pip install pytest pytest-cov
      - run: pytest --cov=compass --cov-report=xml
      - uses: codecov/codecov-action@v3
        with:
          file: ./coverage.xml
```

## GitLab CI

Create `.gitlab-ci.yml`:

```yaml
stages:
  - quality
  - test

code_quality:
  stage: quality
  image: python:3.10
  before_script:
    - pip install -r requirements-dev.txt
  script:
    - python scripts/run_all_checks.py --ci
  artifacts:
    when: always
    paths:
      - lint_reports/
    expire_in: 1 week

tests:
  stage: test
  image: python:3.10
  before_script:
    - pip install -r requirements.txt
    - pip install pytest pytest-cov
  script:
    - pytest --cov=compass --cov-report=xml
  coverage: '/TOTAL.*\s+(\d+%)$/'
  artifacts:
    reports:
      coverage_report:
        coverage_format: cobertura
        path: coverage.xml
```

## Jenkins

### Pipeline Script

```groovy
pipeline {
    agent any
    
    stages {
        stage('Code Quality') {
            steps {
                sh '''
                    pip install -r requirements-dev.txt
                    python scripts/run_all_checks.py --ci
                '''
            }
            post {
                always {
                    archiveArtifacts 'lint_reports/**/*'
                    publishHTML([
                        reportDir: 'lint_reports',
                        reportFiles: '*.html',
                        reportName: 'Code Quality Report'
                    ])
                }
            }
        }
        
        stage('Tests') {
            steps {
                sh '''
                    pip install -r requirements.txt
                    pip install pytest pytest-cov
                    pytest --cov=compass --cov-report=xml
                '''
            }
            post {
                always {
                    publishCoverage(
                        adapters: [
                            coberturaAdapter('coverage.xml')
                        ],
                        sourceFileResolver: sourceFiles('STORE_LAST_BUILD')
                    )
                }
            }
        }
    }
}
```

## Pre-commit Hooks

### Setup pre-commit

1. Install pre-commit:
   ```bash
   pip install pre-commit
   ```

2. Create `.pre-commit-config.yaml`:
   ```yaml
   repos:
     - repo: https://github.com/psf/black
       rev: 23.0.0
       hooks:
         - id: black
           language_version: python3.10
   
     - repo: https://github.com/pycqa/flake8
       rev: 6.0.0
       hooks:
         - id: flake8
           args: [--max-line-length=100]
   
     - repo: https://github.com/pycqa/bandit
       rev: 1.7.5
       hooks:
         - id: bandit
           args: [-r, -ll]
   
     - repo: local
       hooks:
         - id: mypy
           name: mypy
           entry: mypy
           language: system
           types: [python]
           args: [compass, services]
   ```

3. Install hooks:
   ```bash
   pre-commit install
   ```

### Custom Pre-commit Hook

Create `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Run code quality checks before commit

python scripts/run_all_checks.py --ci --no-tests

if [ $? -ne 0 ]; then
    echo "Code quality checks failed. Please fix issues before committing."
    exit 1
fi
```

Make it executable:
```bash
chmod +x .git/hooks/pre-commit
```

## Azure Pipelines

Create `azure-pipelines.yml`:

```yaml
trigger:
  branches:
    include:
      - main
      - develop

pool:
  vmImage: 'ubuntu-latest'

steps:
  - task: UsePythonVersion@0
    inputs:
      versionSpec: '3.10'
    displayName: 'Use Python 3.10'

  - script: |
      pip install -r requirements-dev.txt
    displayName: 'Install dependencies'

  - script: |
      python scripts/run_all_checks.py --ci
    displayName: 'Run code quality checks'
    continueOnError: true

  - task: PublishTestResults@2
    condition: always()
    inputs:
      testResultsFiles: 'lint_reports/**/*.xml'
      testRunTitle: 'Code Quality Reports'

  - task: PublishPipelineArtifact@1
    condition: always()
    inputs:
      targetPath: 'lint_reports'
      artifact: 'lint-reports'
```

## CircleCI

Create `.circleci/config.yml`:

```yaml
version: 2.1

jobs:
  quality:
    docker:
      - image: cimg/python:3.10
    steps:
      - checkout
      - run:
          name: Install dependencies
          command: |
            pip install -r requirements-dev.txt
      - run:
          name: Run quality checks
          command: python scripts/run_all_checks.py --ci
      - store_artifacts:
          path: lint_reports

workflows:
  version: 2
  quality_and_test:
    jobs:
      - quality
```

## Best Practices

1. **Run checks on every PR**: Ensure code quality before merge
2. **Fail fast**: Use `--ci` mode to fail immediately on errors
3. **Store artifacts**: Save reports for later review
4. **Separate concerns**: Run formatting, style, and tests separately if needed
5. **Cache dependencies**: Speed up CI runs by caching pip packages
6. **Parallel execution**: Run different checks in parallel when possible

## Troubleshooting

### Exit codes not working

Ensure you're using `--ci` flag:
```bash
python scripts/run_all_checks.py --ci
```

### Reports not generated

Check that `lint_reports/` directory exists and is writable:
```bash
mkdir -p lint_reports
chmod 755 lint_reports
```

### Missing tools in CI

Install all dev dependencies:
```bash
pip install -r requirements-dev.txt
```

## Examples

See the following files for examples:
- `scripts/run_all_checks.py` - Main check script with CI support
- `scripts/run_lint.bat` - Windows batch script wrapper
- `.github/workflows/` - GitHub Actions workflows (if created)











