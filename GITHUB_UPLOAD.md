# GitHub upload instructions

Recommended repository name:

```text
crbaek/nonconvex-isotonic
```

From the repository root:

```bash
git init
git add README.md CITATION.cff paper_citation.bib MANUSCRIPT_CODE_MAP.md GITHUB_UPLOAD.md LICENSE_MIT_TEMPLATE.txt .gitignore R scripts data results
git commit -m "Initial simulation replication materials"
git branch -M main
git remote add origin git@github.com:crbaek/nonconvex-isotonic.git
git push -u origin main
```

If using HTTPS rather than SSH:

```bash
git remote add origin https://github.com/crbaek/nonconvex-isotonic.git
```

Before making the repository public, confirm the license choice.  If MIT is acceptable, rename `LICENSE_MIT_TEMPLATE.txt` to `LICENSE`.
