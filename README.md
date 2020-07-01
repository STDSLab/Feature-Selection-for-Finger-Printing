# Feature-Selection-for-Finger-Printing
This repository provides MATLAB implementations for our various feature selection methods for FC fingerprinting. This also includes our MATLAB implentation for FC fingerprinting using the selected features. The feature selection methods available in this repository are:

## Node Scoring

### krNodeLearnACSC
Scores node features using Average Cross-Session Cost (ACSC). The ACSC is the average difference in distance between sessions. We compute this by taking the difference between the off-diagonal values and the diagonal-values in the scan session distance matrix. This scan session distance matrix is computed by computing the distance between each FC in set 1 and each FC in set 2.

### krNodeLearnRS
Scores node features using Rank Sum (RS) Cost. The RS cost is computed by summing distance ranks between scans from the same subject. 

### krNodeAnalysis
Using node feature scores, selects the best nodes and performs FC fingerprinting using only the selected node edges.

## Edge Scoring

### edgeLearnACSC
Scores edge features using ACSC.

### edgeLearnRS
Scores edge features using RS cost.

### edgeAnalysis
Using edge feature scores, selects the best edges and performs FC fingerprinting using only the selected edges.

## S Matrix (cluster) Edge Scoring

### mtfLearnACSC
Scores cluster edge features using ACSC.

### mtfLearnRS
Scores cluster edge features using RS cost.

### mtfAnalysis
Using cluster edge feature scores, selects the best cluster edges and performs FC fingerprinting using only the selected cluster edges.
