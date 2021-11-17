## Project description:
This is a registration algorithm by landmark selection and modification.

## Project components:
The algorithm includes three main steps:
1. landmark extraction and matching: Nuclei or corner points are used for landmarks. Landmarks are matched by normalized cross-correlation and partitioned into strong and weak by distance constraints.
2. spatial shift modification: A non-linear weight is used to modify the shift of weak landmarks based on strong matched landmarks.
3. selection by texture and spatial proximity: Modified weak landmarks are selected by distance constraint and SURF features distance.

## Evaluation:
The algorithm provides following metrics to evaluate the registration quality:
- Correlation (COR)
- Mutual Information (MI)
- Mean Squared Error (MSE).

The toolbox structure is given as follows:

    - registration algorithm       ./src/register.m
    - MI calculation               ./MutualInformation.m
    - example                      ./example.m
   