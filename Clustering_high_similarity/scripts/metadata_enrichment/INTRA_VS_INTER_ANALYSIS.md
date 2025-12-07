# Understanding Intra-Component vs Inter-Component Analysis

## Overview

The toolkit performs **TWO complementary types of analysis** that answer different questions about your component metadata patterns.

---

## Intra-Component Analysis (Homogeneity)

### Question
**"Are samples WITHIN each component similar to each other?"**

### What It Measures
- Homogeneity WITHIN components
- Whether samples in the same component share metadata values
- Example: Do all samples in Component 0 come from the same sequencing center?

### Implementation
Done in `analyze_component_metadata.py`:
- Calculates entropy for each field within each component
- Measures frequency of most common value
- Identifies outliers (samples different from component majority)

### Example Results
```
Component 0 (size=150):
  center_name: 98% "BGI" (highly homogeneous)
  bioproject: 100% "PRJNA123" (perfectly homogeneous)
  country: 95% "USA" (highly homogeneous)
  
Component 1 (size=200):
  center_name: 97% "EBI" (highly homogeneous)
  bioproject: 100% "PRJNA456" (perfectly homogeneous)
  country: 92% "China" (highly homogeneous)
```

### Interpretation
- **High homogeneity = Expected for duplicates!**
- Samples cluster because they share metadata (technical replicates, same study, etc.)
- Outliers are interesting (e.g., 5% from different country)

---

## Inter-Component Analysis (Enrichment)

### Question
**"Do DIFFERENT components have DIFFERENT metadata profiles?"**

### What It Measures
- Variation BETWEEN components
- Whether certain metadata values are associated with specific components
- Example: Is Component 0 enriched for Center A while Component 1 is enriched for Center B?

### Implementation
Done in `statistical_enrichment_analysis.py`:
- Creates contingency tables: components × metadata values
- Performs chi-square tests
- Calculates effect sizes (Cramér's V)
- Tests if component membership predicts metadata

### Example Results

**Contingency Table (center_name):**
```
               BGI    EBI    JGI
Component_0    147     2      1     ← Enriched for BGI
Component_1      3   194      3     ← Enriched for EBI  
Component_2      1     4    195     ← Enriched for JGI
```

**Statistical Test:**
- Chi-square: p < 0.001 (highly significant)
- Cramér's V: 0.89 (very strong effect)
- **Interpretation:** Strong association between components and centers

### When This Happens
Two scenarios can produce the same intra-component homogeneity:

**Scenario A: High Inter-Component Variation** (Chi-square significant)
```
Component 0: 100% from Center A
Component 1: 100% from Center B
Component 2: 100% from Center C
```
- High intra-component homogeneity ✓
- High inter-component variation ✓
- Strong statistical enrichment ✓

**Scenario B: Low Inter-Component Variation** (Chi-square NOT significant)
```
Component 0: 100% from Center A
Component 1: 100% from Center A
Component 2: 100% from Center A
```
- High intra-component homogeneity ✓
- Low inter-component variation ✗
- No statistical enrichment ✗

Both have perfect intra-component homogeneity, but only Scenario A shows inter-component enrichment!

---

## Why Both Are Important

### Intra-Component Analysis Tells You:
- ✅ Components are coherent (samples within components are similar)
- ✅ Which metadata fields explain clustering
- ✅ Where outliers exist (samples that don't fit the pattern)

### Inter-Component Analysis Tells You:
- ✅ Different components represent different subpopulations
- ✅ Which metadata fields discriminate between components
- ✅ Statistical significance of associations

---

## Common Patterns

### Pattern 1: Perfect Segregation
**Intra:** High homogeneity (each component dominated by one value)  
**Inter:** Strong enrichment (different components have different values)  
**Example:** Each component is one sequencing center's samples
```
Component 0: 100% BGI
Component 1: 100% EBI
Component 2: 100% JGI
→ Perfect intra-component homogeneity
→ Perfect inter-component enrichment
```

### Pattern 2: Universal Homogeneity
**Intra:** High homogeneity (each component dominated by one value)  
**Inter:** No enrichment (all components have same value)  
**Example:** All samples are human gut metagenomes
```
Component 0: 100% organism="human gut metagenome"
Component 1: 100% organism="human gut metagenome"
Component 2: 100% organism="human gut metagenome"
→ Perfect intra-component homogeneity
→ NO inter-component enrichment (all same!)
```

### Pattern 3: Mixed Within, Different Between
**Intra:** Low homogeneity (components are diverse internally)  
**Inter:** Moderate enrichment (some association exists)  
**Example:** Components mix countries but in different ratios
```
Component 0: 60% USA, 40% China
Component 1: 40% USA, 60% China
→ Moderate intra-component diversity
→ Weak inter-component enrichment
```

### Pattern 4: Complete Randomness
**Intra:** Low homogeneity (components are diverse internally)  
**Inter:** No enrichment (no association)  
**Example:** Random clustering with no metadata signal
```
Component 0: 33% USA, 33% China, 34% UK
Component 1: 32% USA, 35% China, 33% UK
→ Low intra-component homogeneity
→ NO inter-component enrichment
```

---

## The Sparsity "Problem"

### Why Contingency Tables Get Sparse

When you have:
- 100 components
- 50 possible sequencing centers
- Each component homogeneous for ONE center

You get a 100×50 contingency table with:
- 100 non-zero cells (one per component)
- 4,900 zero cells
- **98% sparsity!**

### Why This Is Actually GOOD

**The sparsity IS the pattern!** It means:
- Each component is dominated by ONE value (intra-homogeneity)
- Different components have different values (inter-variation)
- Strong, clear associations (perfect enrichment)

### Old vs New Threshold

**Old approach (too conservative):**
```python
if sparsity > 0.8:  # 80% zeros
    skip_analysis()
```
- Skips many meaningful patterns
- Confuses signal (homogeneity) with noise (missing data)

**New approach (smarter):**
```python
if sparsity > 0.95 and pct_components_with_data < 0.5:
    skip_analysis()
```
- Only skips truly empty fields (missing data)
- Allows sparse but meaningful patterns (high homogeneity)
- Warns when sparse but continues analysis

---

## Workflow Interpretation

### Step 1: Intra-Component Analysis
```
Running: analyze_component_metadata.py

Results:
  center_name: 95% of components have >0.8 homogeneity
  bioproject: 98% of components have >0.9 homogeneity
  country: 85% of components have >0.7 homogeneity
```

**What this means:**
- ✅ Components are coherent
- ✅ Samples within components are similar
- ✅ Strong metadata signals exist

### Step 2: Inter-Component Analysis
```
Running: statistical_enrichment_analysis.py

Results for center_name:
  Sparse table (98% zeros) but 100% components have data, continuing...
  Chi-square: p < 0.001
  Cramér's V: 0.92 (very strong)
  
Results for country:
  Sparse table (95% zeros) but 90% components have data, continuing...
  Chi-square: p = 0.234 (not significant)
  Cramér's V: 0.15 (weak)
```

**What this means for center_name:**
- ✅ High sparsity due to homogeneity (good!)
- ✅ Strong statistical enrichment
- ✅ Different components have different centers
- **Interpretation:** Components segregate by sequencing center

**What this means for country:**
- ✅ High sparsity due to homogeneity (good!)
- ✗ No statistical enrichment
- ✗ Most components have same country
- **Interpretation:** Components are homogeneous for country, but it doesn't discriminate between components (e.g., all USA)

---

## Summary Table

| Analysis | Question | Example | Output |
|----------|----------|---------|--------|
| **Intra-component** | Within-component similarity | "Are samples in Component 0 similar?" | Entropy, frequency, outliers |
| **Inter-component** | Between-component differences | "Is Component 0 different from Component 1?" | Chi-square, p-value, Cramér's V |

| Metric | Intra-Component | Inter-Component |
|--------|-----------------|-----------------|
| **High homogeneity + High enrichment** | Components coherent | Components segregate by metadata |
| **High homogeneity + No enrichment** | Components coherent | All components share same metadata |
| **Low homogeneity + High enrichment** | Components diverse | But some metadata association exists |
| **Low homogeneity + No enrichment** | Components diverse | Random clustering |

---

## Key Takeaway

**For duplicate/identical samples (Jaccard=1.0):**

You should expect:
1. **High intra-component homogeneity** - samples cluster for a reason
2. **Variable inter-component enrichment** - depends on whether reason differs between components

Example:
- All samples in Component 0 from Study A → high intra-homogeneity ✓
- All samples in Component 1 from Study B → high intra-homogeneity ✓
- Study A ≠ Study B → strong inter-enrichment ✓

vs.

- All samples in Component 0 from human gut → high intra-homogeneity ✓
- All samples in Component 1 from human gut → high intra-homogeneity ✓
- All human gut → no inter-enrichment ✗ (but that's OK! It's still informative)

---

**Version:** 1.3  
**Date:** December 4, 2025  
**Purpose:** Clarify the difference between intra and inter-component analyses
