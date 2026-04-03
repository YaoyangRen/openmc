---
name: Beta_eff Benchmark Analysis
description: "Use when you want a fixed template for comparing beta_eff results against MCNP references, analyzing method A versus method B, ranking likely physical causes, and returning a standardized diagnosis report. Keywords: beta_eff, MCNP, benchmark, method A, method B, discrepancy, physical modeling diagnosis"
argument-hint: "Fill in the OpenMC beta_eff result, MCNP reference, error size, files, and any specific concern about method A or method B."
agent: "MCNP Error Analysis"
---
Use this prompt to perform a standardized beta_eff benchmark analysis against MCNP.

Follow the input template exactly. If some fields are unknown, keep the label and write "unknown".

## Input Template

### Case Summary
- Metric: beta_eff
- OpenMC result:
- MCNP reference:
- Error definition: relative error | absolute error | both
- Observed error magnitude:

### Method Comparison
- Method A result:
- Method B result:
- Which method is currently used as the main result:
- What needs explanation: why method B is closer | why method B is farther | compare both methods | unknown

### Inputs and Evidence
- Relevant source files:
- Relevant documents:
- Relevant output files or HDF5 files:
- Known benchmark conditions:
- Known delayed-group convention: 6-group | 8-group | unknown
- Known adjoint quantity meaning: true adjoint flux | scalar importance | unknown

### User Focus
- Priority: physical modeling diagnosis
- Exclude code modification: yes
- Need first function to change later: yes
- Extra concern to check first:

## Required Analysis Rules
- Prioritize physical modeling diagnosis over generic bug hunting.
- Explain the beta_effective computation chain in src before ranking causes.
- Explicitly compare method A and method B against the MCNP reference.
- State whether the better agreement of method B is physically meaningful or likely due to compensating approximations.
- Rank causes by evidence strength, not by brainstorming breadth.
- If evidence is missing, say exactly what is missing.
- Do not modify code.

## Output Format

### Problem
- Metric being compared
- OpenMC result
- MCNP reference
- Error definition
- Observed error magnitude

### Computation Chain
- Main entry point
- Main numerator and denominator functions
- Group metadata or mapping functions involved
- Which implementation branch feeds the reported main result

### Method A vs Method B
- Which method is closer to MCNP
- Why that difference is physically plausible or suspicious
- Whether method B looks like a better model or a compensating approximation

### Findings
- List findings from most likely to least likely
- For each finding include: cause hypothesis, supporting evidence, affected files or functions, confidence

### First Function to Change
- Name one function only
- Explain why it is the best first leverage point

### Modification Plan
- High-priority physical-model changes
- Secondary validation checks
- Risks or assumptions behind the plan

### Missing Information
- List only the missing evidence needed to reduce uncertainty

### Verdict
- One-sentence statement of the most probable primary cause

## Ready-to-Fill Example

### Case Summary
- Metric: beta_eff
- OpenMC result: 0.00584
- MCNP reference: 0.00621
- Error definition: relative error
- Observed error magnitude: -5.96%

### Method Comparison
- Method A result: 0.00531
- Method B result: 0.00584
- Which method is currently used as the main result: method B
- What needs explanation: why method B is closer

### Inputs and Evidence
- Relevant source files: src/beta_effective.cpp, include/openmc/beta_effective.h
- Relevant documents: markdown/BETA_EFFECTIVE_GUIDE.md
- Relevant output files or HDF5 files: beta_eff.h5, flux_mesh.h5, adjoint_flux.h5
- Known benchmark conditions: unknown
- Known delayed-group convention: 6-group
- Known adjoint quantity meaning: scalar importance

### User Focus
- Priority: physical modeling diagnosis
- Exclude code modification: yes
- Need first function to change later: yes
- Extra concern to check first: whether delayed-group mapping and denominator normalization are consistent with MCNP