# Notes on MLEs for the manuscript

This note is for internal use. It records what is mathematically relevant about the MLEs, and also your editorial decisions about what should and should not be said later in the AoS text.

## 1. Editorial decisions to preserve

These are the current decisions and should be treated as binding unless you later change them.

- Keep the MLE discussion short in the paper.
- Do not emphasize multistart.
- Do not emphasize optimizer details such as BFGS or L-BFGS-B.
- For logistic Gaussian, when the real-data analysis is written, simply state the MLE formula.
- For vMF, do not add anything here, because it is already discussed earlier in the paper.
- For HvMF, if it is mentioned, give the formula directly and cite Jensen.
- For JP, keep the delicate computational issues as an internal note for the directors, not as paper text.
- For spherical Cauchy, explain why using the product parameter is mathematically legitimate if we ever mention it.
- For normal and vMF, do not say anything in the new MLE discussion.
- Avoid developing the weighted-MLE versions in the paper unless they are really needed.
- For the small-circle mixtures, if there is nothing beyond basic vlikelihood optimization that matters mathematically, do not say anything.
- For cardioid, if it is just likelihood maximization, do not say anything.
- Ignore the rotational logit-normal mixture for now.
- Leave the axial truncated-normal mixture for later.

## 2. Global mathematical message

The clean global statement is the following.

- In every composite scenario, the unknown parameters are estimated by maximum likelihood.
- In the models with a closed-form characterization, that exact characterization is used.
- In the models without a closed form, the estimator is the joint numerical maximizer of the log-likelihood.
- When a computation can be written in successive steps, this is only because the joint likelihood equations reduce to that form; it is not a different two-stage estimator.

So, on the specific point you raised: the important distinction is not "joint" versus "successive" at the level of implementation, but whether the successive expression still equals the exact joint MLE. In the standard closed-form cases here, the answer is yes.

## 3. Joint versus successive

There are two different notions.

- First, one may solve first for one block because the likelihood equations force that reduction. This is still the joint MLE.
- Second, one may estimate one block and then another from a different criterion. That would be a genuinely sequential estimator.

What we have in the relevant closed-form examples is the first situation, not the second.

## 4. Model-by-model summary

| Model | Unknown parameters | Joint or successive? | How is the MLE obtained? | What is relevant to say? |
|---|---|---|---|---|
| Normal on $\mathbb{R}$ | $\mu$, $\sigma$, or both | Joint. It can be written in successive form, but it is still the joint MLE. | If only $\mu$ is unknown, then $\hat{\mu}=\frac{1}{n}\sum_{i=1}^n X_i$. If only $\sigma$ is unknown with $\mu$ fixed, then $\hat{\sigma}^2=\frac{1}{n}\sum_{i=1}^n (X_i-\mu)^2$. If both are unknown, then $\hat{\mu}=\frac{1}{n}\sum_{i=1}^n X_i$ and $\hat{\sigma}^2=\frac{1}{n}\sum_{i=1}^n (X_i-\hat{\mu})^2$. | Closed form. Since you do not want to discuss normal here, this is only background. |
| Logistic Gaussian on the simplex | $\mu_{\mathrm{ilr}}$, $\Sigma_{\mathrm{ilr}}$, or both | Joint after the ilr transform. The covariance uses $\hat{\mu}_{\mathrm{ilr}}$, but this is still the exact joint Gaussian MLE in ilr coordinates. | Let $Z_i=\operatorname{ilr}(X_i)$. Then $\hat{\mu}_{\mathrm{ilr}}=\frac{1}{n}\sum_{i=1}^n Z_i$ and $\hat{\Sigma}_{\mathrm{ilr}}=\frac{1}{n}\sum_{i=1}^n (Z_i-\hat{\mu}_{\mathrm{ilr}})(Z_i-\hat{\mu}_{\mathrm{ilr}})^{\top}$. <span style="color:red">Weighted versions can be written analogously with weighted averages if you later decide that the bootstrap version should be made explicit.</span> | Closed form after an isometric transform to Euclidean space. If mentioned, the right wording is that this is just the Gaussian MLE in ilr coordinates. |
| vMF on $S^q$ | canonical parameter $\xi=\kappa\mu$ | Joint. The decomposition into direction and concentration is equivalent to the joint MLE. | The score equation reduces to the sample resultant. | Already treated elsewhere in the paper. Do not add a new MLE discussion here. |
| HvMF on $H^2$ | $\xi$ and $\kappa$ | Joint, with exact closed form. | If $S=\sum_{i=1}^n X_i$ and $R=\sqrt{-\langle S,S\rangle_M}$, the MLE exists when $R>n$ and is given by $\hat{\xi}=S/R$ and $\hat{\kappa}=n/(R-n)$. <span style="color:red">In the weighted implementation, with $S=\sum_i w_i X_i$ and $W=\sum_i w_i$, the same formula becomes $\hat{\xi}=S/R$ and $\hat{\kappa}=W/(R-W)$.</span> | This is the main additional formula worth recording, provided it is cited directly to Jensen. |
| Spherical Cauchy on $S^2$ | $\mu$ and $\rho$ | Joint numerical MLE. | The likelihood is rewritten in terms of $\phi=\rho\mu$, with $\|\phi\|<1$, and then optimized numerically after an unconstrained reparametrization. | If ever mentioned, the mathematically interesting point is not only legitimacy but utility: the pair $(\mu,\rho)$ lives on a constrained space, while $\phi$ lives in the open unit ball and can then be mapped smoothly to an unconstrained Euclidean variable. |
| Jones-Pewsey on S^2 | mu, kappa, psi | Joint numerical MLE. | There is no closed-form formula in the implemented parametrization; the estimator is obtained by numerical maximization of the joint log-likelihood. | Keep the computational complications out of the paper. Use them only as an internal note if needed. |
| Small circle on S^2 | mu, kappa, nu | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | Unless later needed for a specific reason, say nothing. |
| Symmetric two-small-circle mixture on S^2 | mixture parameters | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | Unless later needed for a specific reason, say nothing. |
| Weighted two-small-circle mixture on S^2 | mixture parameters | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | Unless later needed for a specific reason, say nothing. |
| Cardioid / spherical cardioid | mu and rho | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | If that is all, say nothing. |
| Rotational beta mixture on S^2 | mixture parameters | Joint numerical MLE. | The final estimator is obtained by numerical maximization of the joint log-likelihood. The sample splitting and moment matching only serve to generate starting values. | If it appears in the paper, the honest mathematical description is still "joint numerical MLE". The moment-matching step is not a different estimator. |
| Rotational logit-normal mixture on S^2 | mixture parameters | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | Ignore for now. |
| Axial truncated-normal mixture | mixture parameters in the axial model | Joint numerical MLE. | Numerical maximization of the joint log-likelihood. | Leave for later. |

Except where explicitly marked in red, the formulas above are written in the ordinary unweighted form, because that is the cleaner version to keep in mind for exposition.

## 5. Which models are really different from "just maximizing a log-likelihood"?

If by that phrase one means "the only thing happening is a generic black-box optim call", then the answer is:

- Logistic Gaussian: different. Closed form after ilr transformation.
- HvMF on $H^2$: different. Exact closed form.

The remaining models are indeed numerical maximum-likelihood fits. For paper purposes, the only computational distinctions that seem worth preserving internally are:

- Spherical Cauchy: the unconstrained fit is based on the reparametrization $\phi=\rho\mu$, which replaces the constrained parameter pair $(\mu,\rho)$ by a single Euclidean vector in the open unit ball.
- JP: there are numerical safeguards in the implementation, but this should stay out of the paper text unless it becomes unavoidable.
- Rotational beta mixture: the split-and-match step is only for initialization, not part of the estimator itself.

So, if the question is "is there any MLE worth explicit mention?", the answer is:

- In the current state of the paper, the only new one really worth mentioning is HvMF, if that model is going to be mentioned at all.
- For logistic Gaussian, a direct formula is enough when the real-data analysis is written.
- For vMF, nothing new should be added here.
- For JP, cardioid and the small-circle family, the default should be to say nothing unless some theoretical necessity appears.

## 6. What I would recommend stating explicitly in the paper

The cleanest statement is:

- In every composite-null experiment, all unknown parameters are estimated by maximum likelihood.
- Whenever the joint MLE admits a closed-form characterization, that characterization is used.
- Otherwise, the joint log-likelihood is maximized numerically.

If you want a model-specific version adapted to your current editorial criteria, the safest one is:

- In the logistic Gaussian model, after the ilr transformation, the MLEs are the usual Gaussian estimators for the mean and covariance.
- If HvMF is mentioned, cite Jensen and state the closed-form estimator directly.

## 7. Internal note on JP

This section is not for the paper text.

There are two numerical safeguards in the JP implementation that are worth remembering internally.

- Near the vMF regime, when $|\kappa\psi|$ is very small, the code switches to the exact vMF likelihood for stability. The default threshold is $10^{-3}$.
- In the benchmark and calibration workflows relevant here, the code often enforces the cap $|\kappa\psi|\le 6$.

<span style="color:blue">Note for the directors. Near $\psi=0$, we switch to the exact vMF likelihood when $|\kappa\psi|$ is very small. This is not a change of model, but a numerical safeguard used to avoid computational instability in the JP likelihood near the vMF regime.</span>

Your recollection was correct: the quantity $|\kappa\psi|\le 6$ does appear in the implementation.

## 8. Internal note on spherical Cauchy

If we ever need to justify the parameter $\phi=\rho\mu$, the reason is simple.

- The original parameter is a pair $(\mu,\rho)$ with $\|\mu\|=1$ and $0\le \rho<1$.
- The product $\phi=\rho\mu$ belongs to the open unit ball.
- For $\rho>0$, the correspondence between $(\mu,\rho)$ and $\phi$ is one-to-one, since $\rho=\|\phi\|$ and $\mu=\phi/\|\phi\|$.
- At $\rho=0$, the direction $\mu$ is intrinsically non-identifiable anyway, so nothing is lost by using $\phi$.

Therefore this is not a trick that changes the estimator; it is just a legitimate reparametrization of the same model.

The practical advantage is that it removes the geometric constraints from the optimization.

- In the original variables, one must optimize over a direction $\mu\in S^2$ together with a radial parameter $\rho\in[0,1)$.
- In the variable $\phi$, both parameters are merged into a single vector with $\|\phi\|<1$.
- That open unit ball can then be mapped to an unconstrained $u\in\mathbb{R}^3$ by
  $\phi=u/\sqrt{1+\|u\|^2}$.

So the point is not conceptual elegance; the point is that numerical optimization becomes an ordinary unconstrained Euclidean problem, while still representing exactly the same statistical model.

Also, the analytic gradient was in fact used in the spherical Cauchy workflows, so mentioning that would be true if ever needed.

## 9. Source pointers in the code

The main code locations are:

- Normal: `bootstrap/model_specs.R`, lines 262-321.
- Logistic Gaussian: `bootstrap/model_specs.R`, lines 755-835.
- vMF: `bootstrap/model_specs.R`, lines 1031-1108, and `utils.R`, lines 1478-1486.
- HvMF: `bootstrap/model_specs.R`, lines 1454-1484, and `utils.R`, lines 305-372.
- Spherical Cauchy: `bootstrap/model_specs.R`, lines 1568-1592, and `utils.R`, lines 2216-2352.
- JP: `bootstrap/model_specs.R`, lines 1280-1312, and `utils.R`, lines 2928-3448.
- Small circle: `utils.R`, lines 4966-5144.
- Small-circle weighted mixture: `utils.R`, lines 6186-6270.
- Small-circle symmetric mixture: `utils.R`, lines 6503-6713.
- Beta mixture: `utils.R`, lines 7578-7660.
- Logit-normal mixture: `utils.R`, lines 8519-8602.
- Cardioid: `bootstrap/cardioid_model_spec.R`, lines 191-321.
