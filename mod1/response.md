We thank the referee for providing constructive and thoughtful
comments to improve the manuscript. We have addressed all the concerns
raised by the referee and details can be found in the revised version
-- please see attached. A brief description on how we incorporated the
changes can be found below.


# Major comments
1:
- There is already a lensing model in H14 (a paper the authors refer
  to numerous times) and it is unclear why they redid the lensing
  model (they use the same SMA data as H14 but a different code). If
  the derived model is significantly different from H14 then this
  should be highlighted otherwise the lens modeling section isn't
  necessary in this paper.

> Most parameters obtained in our new analysis with our new lensing
  code are indeed consistent with the earlier model reported by H14,
  such that the results are only slightly improved. However, H14 did
  not report the full range of best-fit parameters and uncertainties
  in their analysis, such that our results are significantly more
  detailed. To take this comment into account, we thus significantly
  shortened the text in this section by removing details on the lens
  modeling approach (and including a reference to the literature
  instead), and we have moved the lens mass and consistent results
  into Table 1 - which is used to briefly summarize the full range of
  parameters. This includes the important inferred half-light radius,
  which is used later in the analysis.

2:
- The SMG is thought to harbor a luminous AGN. In that case, the LIR
  will certainly have a contribution from the AGN that should be
  accounted for otherwise the SFR will be overestimated: see
  e.g. Coppin+2010, Kirkpatrick+2012, Magdis+2013,
  Kirkpatrick+2014. Please correct LIR for the AGN emission prior to
  estimating SFR and include an AGN component in the SED fitting.

> In our analysis, we fit SED models to the rest-frame FIR, which has
  been found to have negligible AGN contribution to rest-frame FIR (≥
  60µm; Carilli+05, Netzer+07, Mullaney+11, Harrison+15). Since we do
  not have MIR measurements solely from the background SMG and due to
  the positive K correction blue-ward of dust peak, the foreground radio
  galaxy contributes a non-negligible amount to the MIR luminosity, we
  therefore didn't fit an AGN component. Including an additional
  component, and thus free parameters to the SED model will be
  over-fitted due to the limited amount of data at hand. Instead our
  SED model includes a power-law to the blue side of the dust peak,
  which accounts for the MIR excess and allows us to estimate the IR
  luminosity (Kirkpatrick+15, Riechers+13, Johnson+13, Casey+12). We
  derived the SFR and other corresponding properties using the L_FIR,
  (45-122.5), which is mostly due to SF-driven dust emission.

3:
- section 4.3 is way too wordy and can be cut down to a few paragraphs
  and a table summarizing these key values (already given in Table
  4). These are all standard calculations and you can refer to CW14
  for the equations rather than repeating them here.

> We removed standard equations from the manuscript, and combined the
  sections reporting SFE, SFR, and SFRSD into one. We also shortened
  the section.

4:
- You should also summarize the evidence for the AGN in this SMG
  (maybe when you mention this in the intro) so the reader doesn't
  have to read H14.

> We rephrased the introduction to show more clearly the evidence of
  an AGN.

5:
- [a further concern/question about Table 4 is whether all the same
  assumptions such as excitation, alpha_co, etc. were made for the
  other SMGs otherwise the comparisons are meaningless]

> We added notes in the caption of Table 4 to recap all assumptions
  used in deriving these quantities.

6:
- Section 5: why do the authors cherry-pick these 2 SMGs to make
  comparisons? there are many other sources that could be chosen with
  a range of properties. A comparison to these 2 SMGs is not useful at
  all given the diversity in this population (e.g. the sentence "This
  demonstrates the large ..." is not true when comparing only 3
  objects). I suggest comparing to the general properties observed for
  SMGs (as summarized in CW14 or Casey+14) which will be more useful
  in placing this object in the context of the SMG population.

> We took into account this comment by adding a compilation of SMG and type-2 
> QSOs properties to Table 4. Please see the table notes and comments for
  details. We have therefore included a comparison to these sample-averaged
  values in our conclusion. 
  We initially picked out two well-studied and apparently similarly bright (submillimeter fluxes >400 mJy)gravitationally lensed
  SMGs, to be able to investigate a broader range of physical properties - some of which are not well-constrained for the broader population. Thus, we decided to also keep this comparison in the paper.

## 7
- the discussion/conclusions should focus on what is new from this
  paper and not what is already presented in H14. E.g. the last
  sentence about this SMG transitioning from starburst to QSO
  phase. How is that conclusion made or strengthened with this new
  data?

> We updated the conclusions.



# Minor comments
1:
- the depletion time seems short compared to other SMGs which are
  ~100Myr although maybe it is consistent with other SMGs with AGN in
  them.

> Greve+05 report a range of depletion time τ = 40-100 Myr for SMGs,
  please see also Table 4.


2:
- section headings in 3 are strange. Remove the words "New results"
  since this is obvious and a requirement for publication. Also remove
  subheading or second line in heading 3.2.

> We removed "New results" in the subject headings.

3:
- the left panels on Fig 1 and 2 should show the same region (have the
  same RA/DEC axes). It is also tough to tell by eye if the VLA9GHz
  and SMA emission is coincident. Perhaps these can be overlaid on
  Fig 1.

> We adjusted the panel sizes as suggested, and marked the locations
  of the peak emission from the SMA as white crosses on the right
  panel in Figure 1.

4:
- figure 3: orange is difficult to see. What is the SNR of the
  residuals in the right panel?

> We modified our caption to state that the quoted contour levels are
  the same for both panels in the figure. We changed the critical line
  color to green for better readability.

5:
- figure 4: again colors (red, orange) are hard to distinguish on a
  print out version of the paper.

> We did not use orange in this figure. Please kindly clarify.

6:
- specify the wavelengths range assumed for FIR.

> We mentioned the FIR range (42.5-122.5 micron) in Section 4.2.2.

7:
- last para section 4.2.2: it is limited data in the submm not the FIR
  [need to be probing Rayleigh Jeans] that makes the dust mass weakly
  constrained. I don't understand what the last sentence in this para
  means - I think there are words missing.

> We clarified this issue by adding the phrase "rest-frame" to FIR. We
  also restructured the last paragraph to convey our message more
  clearly.

8: 
- SF surface density: from resolved studies of SMGs is it reasonable
  to divide half the SFR and gas mass by the area out to the
  half-light radii? give appropriate references to justify this
  assumption.

> We now provide a few citations to show that this is a common
  assumption.


9:
- intro: the first para makes it sound like our knowledge of SMGs
  comes primarily from SMGs found in large sky Herschel surveys and
  SPT but this is not true. Please include references to earlier work
  that established the basic properties of SMGs (redshift
  distribution, luminosities, SFRs, Mgas, etc) prior to Herschel even
  launching.

> We restructured the paragraph and added earlier discoveries.


10:
- the fact that the lensing system has two components isn't presented
  clearly in this paper. I had to read H14 to understand what the
  authors meant when they mentioned companion B. Please make this
  clearer so the reader does not have to read H14

> We added a sentence to the introduction section to elaborate on this, please see the revised manuscript.

11:
- last few paras of paper: what does "breaking up its dust cocoon"
  mean? this is the first time this is mentioned and should be
  explained if it is relevant.

> We mentioned this as a reference from Haas+14. We have removed this
  phrase in the revised version given that this has been mentioned in
  Haas+14, and is not necessary in our discussion.

12:
- end of section 4.3.1: the factor of 2 for dust mass is overshadowed
  by the factor of >5 in alpha_co

> We removed that sentence.

13:
- related to depletion time, there is a new Nature paper by Narayanan
  et al (published Sept 23) which discusses how SMGs in sims can last
  for a Gyr with replenishment from the IGM. In that case, is the
  depletion timescale really useful? what do we learn if the galaxy is
  getting new gas?

>The derived depletion timescale is derived the based on the
 definition of the quantity τ = Mgas/SFR, which describes the time it
 would take for a galaxy to use up the available molecular gas
 reservoir given the current star formation rate, assuming no
 replenishment of gas (Carilli & Walter 2013). Hence, this quantity
 provides a lower limit on the timescale for SF activity. Based on our
 observations, we find no evidence of gas inflow in SMMJ0939.

13:
- dynamical mass calculation: can you estimate how much higher this
  would be given observations showing the different radii between the
  molecular and dust emission in SMGs. How different are the radii?

>The half-light radii of dust-emitting regions in SMGs typically are ~
 1.2 kpc based on a recent study by Simpson+15, who collected ALMA 870
 µm data for 52 SMGs. Tacconi+06 report a median FWHM size of 4.1 ±1.6
 kpc based on high-J CO measurements. Given the small number
 statistics on resolved CO(1-0) emission in QSOs, where Riechers+11
 find no evidence of extent cold gas emission, we did not include an
 additional estimate of the dynamical mass assuming some physical size
 from recent observations (still limited) in the paper. However, we
 have noted these physical sizes in the text, which allows readers to
 scale the dynamical mass correspondingly should they wish to do
 so. Please see the revised manuscript.

14:
- appendix: "no evidence of correlation" do the authors mean
  "covariance" here?

> We used the word "correlation" to indicate that regardless of the
  variation found in dust mass due to a different prior imposed on the
  emissivity, the IR luminosity is insensitive to such changes in
  dust mass. Both the dust mass and IR luminosity are inferred from the
  best-fit SED model, and are not parameters in the model (and the
  loss function). For these reasons, we do not mean the covariance
  here.

