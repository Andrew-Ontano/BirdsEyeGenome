# Agent Instructions & Project Context

> This document defines the operational boundaries, architecture guidelines, and coding standards for autonomous agents working within this repository.

---

## 1. Project Overview & Architecture
* **Core Purpose:** Birds Eye Genome takes in genomic data, outputs short length k-mer counts broken wby windows, then uses those k-mers to visualize sequence structure.
* **Tech Stack:** Python 3.12, JavaScript, Conda, SQLite
* **Directory Structure:**
  * `/src` - Core pipeline logic and modular stages
  * `/tests` - Unit, integration, and end-to-end test suites
  * `/configs` - Environment variables and pipeline configuration schemas

---

## 2. Agent Operational Boundaries & Safety
* **Autonomous Permissions:**
  * Allowed to create new test files, modify implementation files inside `/src`, and run test commands.
  * **Prohibited:** Never alter CI/CD GitHub action workflows without explicit confirmation, or commit directly to `main`.
* **Error Handling Policy:** If a test fails or a dependency conflict arises, pause, analyze the traceback, and propose a fix rather than attempting destructive overrides (e.g., forcing package installations).

---

## 3. Code Quality & Style Standards
* **Language/Framework Conventions:**
  * Strict typing required (e.g., Python type hints for all function signatures. OR operands should not be used, accept the singular type hint as stable).
  * Use descriptive variable names; avoid abbreviations unless widely accepted in domain context.
* **Documentation:** Every public module, class, and core pipeline function must include docstrings explaining parameters, return types, and side effects.
* **Formatting & Linting:** Code must comply with project linters (e.g., `ruff`, `black`, or `flake8`). Always run the linter before completing a task.

---

## 4. Testing & Validation Requirements
* **Test-Driven Mindset:** Whenever a new feature or pipeline stage is added, corresponding unit tests must be written or updated in `/tests`.
* **Execution Commands:** Agents must use the following commands to validate changes:
  * Run unit tests: `pytest tests/unit`
  * Run lint check: `ruff check .`

---

## 5. Workflow Execution Steps (For Jules Pipelines)
When assigned a pipeline task or issue, follow this sequence:
1. **Analyze:** Review the target requirements, existing code patterns, and relevant test files.
2. **Implement:** Write clean, modular code conforming to the style guidelines above.
3. **Verify:** Run the designated test and lint commands to ensure zero regressions.
4. **Report:** Provide a concise summary of changes made, files touched, and validation results.