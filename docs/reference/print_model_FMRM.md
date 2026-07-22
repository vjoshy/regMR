# Printing Helper Function for FMRM()

Printing Helper Function for FMRM()

## Usage

``` r
print_model_FMRM(
  selected_parameters,
  selected_compartment,
  family,
  information_criteria
)
```

## Arguments

- selected_parameters:

  parameters for finite mixture regression model, depends on family.

- selected_compartment:

  An integer greater than or equal to one representing the number of
  mixture components (groups) in a finite mixture regression model.

- family:

  A string of characters specifying the distribution of the finite
  mixture regression model being fit to the data. Parameter updates are
  altered depending on the inputted family.

- information_criteria:

  A string of characters specifying the information criteria for model
  selection purposes. The model that minimizes the information criteria
  over all group counts and lambda-alpha pairs will be selected. Current
  accepted types include the default Bayesian Information Criterion
  (BIC) ("bic"), group-structured Extended BIC (gEBIC) ("gebic"), Akaike
  Information Criterion (AIC) ("aic"), and Integrated Classification
  Likelihood (ICL) Criterion ("icl").

## Value

No return value, called for side effects
