---
layout: page
redirect_from: 
---

# TCP files

| Files       |  Last update |
|:-------------|:------------------|
| <a href="./TCP_Supplementary_Appendix.pdf" target="_blank">TCP_Supplementary_Appendix.pdf</a>           | 30/01/2019 |
| <a href="./TCP_Paper.pdf" target="_blank">TCP_Paper.pdf</a> | 30/01/2019 |


##### Abstract 
What is the relationship between international trade and business cycle synchronization? Using data from 40 countries, we find that GDP comovement is significantly associated with trade in intermediate inputs but not with trade in final goods. Motivated by this new fact, we build a model of international trade that is able to replicate the empirical trade-comovement slope, offering the first quantitative solution for the *Trade Comovement Puzzle*. The model relies on (i) global value chains, (ii) price distortions due to monopolistic competition and (iii) fluctuations in the mass of firms serving each country. The combination of these ingredients creates a link between domestic measured productivity and foreign shocks through trade linkages, generating a disconnect between technology and measured productivity. Finally, we provide empirical evidence for the importance of these elements in generating a link between foreign shocks and domestic GDP.

##### TCP estimates
<span id="sensitive_TC_slope"
label="sensitive_TC_slope">\[sensitive\_TC\_slope\]</span>

<table>
<caption>Sensitive analysis: TC-slope</caption>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: left;">Coefficient on trade in inputs</th>
<th style="text-align: left;">Coefficient on trade in final goods</th>
<th style="text-align: left;">GDP Filter</th>
<th style="text-align: left;">Countries | Obs.</th>
<th style="text-align: center;">Period</th>
<th style="text-align: center;">TW</th>
<th style="text-align: center;">CP</th>
<th style="text-align: center;"></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Whole Sample</td>
<td style="text-align: left;">0.053<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.030</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 2,900</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">20 years TW</td>
<td style="text-align: left;">0.074<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.054</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 1,450</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Excluding EU CP</td>
<td style="text-align: left;">0.056<span class="math inline">**</span></td>
<td style="text-align: left;">0.005</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 2,280</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">Excluding USSR</td>
<td style="text-align: left;">0.064<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.006</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">34 | 2,244</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">Alternative TW</td>
<td style="text-align: left;">0.081<span class="math inline">***</span></td>
<td style="text-align: left;">0.014</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">34 | 2,244</td>
<td style="text-align: center;">1970-1999</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">4Digits SITC</td>
<td style="text-align: left;">0.058<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.045<span class="math inline">*</span></td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">36 | 2,520</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">ISIC classification</td>
<td style="text-align: left;">0.059<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.045<span class="math inline">*</span></td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">36 | 2,520</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;">1Digit Agg. sectors</td>
<td style="text-align: left;">0.088</td>
<td style="text-align: left;"><span class="math inline">-</span>0.044</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">38 | 1,291</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: left;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;"><span class="math inline">level(trade)</span></td>
<td style="text-align: left;">33,96<span class="math inline">*</span></td>
<td style="text-align: left;"><span class="math inline">-</span>34.92</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 2,900</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;"><span class="math inline">log(mean(trade))</span></td>
<td style="text-align: left;">0.044<span class="math inline">*</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.027</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 2,900</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="even">
<td style="text-align: left;"><span class="math inline">max{T/GDP_1,T/GDP_2}</span></td>
<td style="text-align: left;">0.052<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.032</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">40 | 2,900</td>
<td style="text-align: center;">1970-2009</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
<tr class="odd">
<td style="text-align: left;">STAN data</td>
<td style="text-align: left;">0.209<span class="math inline">**</span></td>
<td style="text-align: left;"><span class="math inline">-</span>0.107</td>
<td style="text-align: left;">HP</td>
<td style="text-align: left;">20 | 760</td>
<td style="text-align: center;">1995-2014</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;">Yes</td>
<td style="text-align: center;"></td>
</tr>
</tbody>
</table>
