We thank the referee for providing constructive and thoughtful comments to improve the manuscipt. We have addressed all the concerns raised by the referee and details can be found in the revised version -- please see attached. A brief description on how we incoprated the changes can be found below. 


# Major comments
1:
- There is already a lensing model in H14 (a paper the authors refer to numerous times) and it is unclear why they redid the lensing model (they use the same SMA data as H14 but a different code). If the derived model is significantly different from H14 then this should be highlighted otherwise the lens modeling section isn't necessary in this paper.
> We moved the lens mass and consistent results into table 1 and its caption, and removed details of the lens model and shortened the section. However, we retained this section in this paper as we provide the best-fit parameters that were not reported in Haas+14 and since we refered to the inferred half-light radius later in the analysis section.

2:
- The SMG is thought to harbor a luminous AGN. In that case, the LIR will certainly have a contribution from the AGN that should be accounted for otherwise the SFR will be overestimated: see e.g. Coppin+2010, Kirkpatrick+2012, Magdis+2013, Kirkpatrick+2014. Please correct LIR for the AGN emission prior to estimating SFR and include an AGN component in the SED fitting. 
> In our analysis, we fit SED models to the rest-frame FIR, which has been found to have negligible AGN contribution to rest-frame FIR (≥ 60µm; Carilli+05, Netzer+07, Mullaney+11, Harrison+15). Since we do not have MIR measurements solely from the background SMG and due to the postive K corretion blueward of dust peak, the forground radio galaxy contributes a non-negligible amount to the MIR luminosity, we therefore didn't fit an AGN component. Including an additional component, and thus free parameters to the SED model will be over-fitted due to the limited amount of data at hand. Instead our SED model includes a power-law to the blue side of the dust peak, which accounts for the MIR excess and allows us to estimate the IR luminosity (Kirkpatrick+15, Riechers+13, Johnson+13, Casey+12). We derived the SFR and other corresponding properties using the L_FIR, (45-122.5), which is mostly due to SF-driven dust emission.

3:
- section 4.3 is way too wordy and can be cut down to a few paragraphs and a table summarizing these key values (already given in Table 4). These are all standard calculations and you can refer to CW14 for the equations rather than repeating them here.
> We removed standard equations from the manuscript, combined the sections reporting SFE and SFR, SFRSD into one, and trimmed a few sentences.

4:
- You should also summarize the evidence for the AGN in this SMG (maybe when you mention this in the intro) so the reader doesn't have to read H14.
> We rephrased the introduction to show more clearly the evidence of an AGN.

5:
- [a further concern/question about Table 4 is whether all the same assumptions such as excitation, alpha_co, etc. were made for the other SMGs otherwise the comparisons are meaningless]
> We added notes in the caption below Table 4 to recap assumptions used in deriving these quantities.

## 6 (NEED TO WORK ON):
- Section 5: why do the authors cherry-pick these 2 SMGs to make comparisons? there are many other sources that could be chosen with a range of properties. A comparison to these 2 SMGs is not useful at all given the diversity in this population (e.g. the sentence "This demonstrates the large ..." is not true when comparing only 3 objects). I suggest comparing to the general properties observed for SMGs (as summarized in CW14 or Casey+14) which will be more useful in placing this object in the context of the SMG population. 

> We picked out two gravitaitonally lensed SMGs to show the ranges of SMGs properties, of which intrinsically sit at faint ends of the luminosity distribution, and thus would not have been studied in sensitivity-limited samples. 
We have added a compilation of SMG and type-2 QSOs properties to the Table, where SMGs properties are derived from Sharon+15 on a sample of SMGs with both CO(3-2) and CO(1-0) line measurements. Hence the SMG properties are free of assumption on the excitation conditions. Please see tablenote and comments for details. We have incoporated a comparison to these values in our conclusion. Please see conclusion section for other changes.

## 7 (NEED TO WORK ON):
- the discussion/conclusions should focus on what is new from this paper and not what is already presented in H14. E.g. the last sentence about this SMG transitioning from starburst to QSO phase. How is that conclusion made or strengthened with this new data?

> We have updated the conclusions.



# Minor comments
1:
- the depletion time seems short compared to other SMGs which are ~100Myr although maybe it is consistent with other SMGs with AGN in them.
> Greve+05 report a range of depletion time τ = 40-100 Myr for SMGs. 


2:
- section headings in 3 are strange. Remove the words "New results" since this is obvious and a requirement for publication. Also remove subheading or second line in heading 3.2.
> We have removed "New results" in the subject headings.

3:
- the left panels on Fig 1 and 2 should show the same region (have the same RA/DEC axes). It is also tough to tell by eye if the VLA9GHz and SMA emission is coincident. Perhaps these can be overlayed on Fig 1.
> We marked the locations of the peak emission from the SMA as white crosses on the right panel on figure 1.

4:
- figure 3: orange is difficult to see. What is the SNR of the residuals in the right panel?
> We modified our caption to state that the quoted contour levels are the same for both panels in the figure. We changed the critical line color to green to easily visualization. 

5:
- figure 4: again colors (red, orange) are hard to distinguish on a print out version of the paper.
> We were unable to make adjustment regarding this as we did not use orange in the figure. Please kindly clarify this. 

6:
- specify the wavelengths range assumed for FIR.
> We mentioned the FIR range in Section 4.2.2.

7:
- last para section 4.2.2: it is limited data in the submm not the FIR [need to be probing Rayleigh Jeans] that makes the dust mass weakly constrained. I don't understand what the last sentence in this para means - I think there are words missing.
> We clarify this by adding the phrase "rest-frame" to FIR. We also restructured the last paragraph to convey our message more clearly.

8: 
- SF surface density: from resolved studies of SMGs is it reasonable to divide half the SFR and gas mass by the area out to the half-light radii? give appropriate references to justify this assumption.
> We provided a few citations to show that this is typical.


9:
- intro: the first para makes it sound like our knowledge of SMGs comes primarily from SMGs found in large sky Herschel surveys and SPT but this is not true. Please include references to earlier work that established the basic properties of SMGs (redshift distribution, luminosities, SFRs, Mgas, etc) prior to Herschel even launching.
> We restructured the paragraph and added earlier discoveries.

10:
- the fact that the lensing system has two components isn't presented clearly in this paper. I had to read H14 to understand what the authors meant when they mentioned companion B. Please make this clearer so the reader does not have to read H14
> We added a sentence to elaborate on this, please see the revised manuscript.

11:
- last few paras of paper: what does "breaking up its dust cocoon" mean? this is the first time this is mentioned and should be explained if it is relevant.
> We mentioned this as a reference from Haas+14. We have removed this phrase in the revised version given that this has been mentioned in Haas+14, and is not necessary in our discussion.

12:
- end of section 4.3.1: the factor of 2 for dust mass is overshadowed by the factor of >5 in alpha_co
> We removed that sentence.

13:
- related to depletion time, there is a new Nature paper by Narayanan et al (published Sept 23) which discusses how SMGs in sims can last for a Gyr with replenishment from the IGM. In that case, is the depletion timescale really useful? what do we learn if the galaxy is getting new gas?
>The derived depletion timescale is derived the based on the definition of the quantity τ = Mgas/SFR, which descibes the time 
it would take for a galaxy to use up the available molecular gas reservoir given the current star formation rate, assuming no 
replenishment of gas (Carilli & Walter 2013). Hence, this quantity provides a lower limit on the timescale for SF activity. 
Based on our observations, we find no evidence of gas inflow in SMMJ0939. 

13:
dynamical mass calculation: can you estimate how much higher this would be given observations showing the different radii between the molecular and dust emission in SMGs. How different are the radii?
>The half-light radius of dust emitting region is ~ 1.2 kpc based on a recent study by Simpson+15, who collected ALMA 870 µm data over 52 SMGs. Tacconi+06 report a median FWHM size of 4.1 ±1.6 kpc based on CO measurements. Furthermore, Ivison+11 found evidence of extended ground state CO emittiong region of ~7 kpc. Please see the revised manuscipt. 
 
**@DR, should I calculate a dynamical mass using Median from Tacconi+06 as well? I get Mdyn ~ 36.8e10 Msun. Should I include these sizes in the paper? or just feed this to the referee.?**

14:
- appendix: "no evidence of correlation" do the authors mean "covariance" here?
> We used the word "correlation" to indicate that regardless of the variation found in dust mass due to a different prior imposed on the emissivity, the IR luminosity is insensivtive to such changes in dust mass. Both the dust mass and IR lumonisty are inferred from the best-fit SED model, and are not parameters in the model (and the loss function). For these reasons, we do not mean the covariance here.










