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
\begin{table}[!htbp] \centering 
\small
  \caption{Sensitive analysis: TC-slope} 
  \label{sensitive_TC_slope} 
\begin{tabular}{@{\extracolsep{0pt}}lp{2.2cm}p{2.2cm}p{1.0cm}p{1.8cm}cccc} 
\toprule & Coefficient on  trade in inputs  & Coefficient on trade in final goods & GDP Filter & Countries | Obs. & Period & TW & CP \\ \midrule 
\multicolumn{8}{l}{\textit{Sample selection}} \\ [0.1cm]
Whole Sample & 0.053$^{**}$ & $-$0.030 & HP & 40 | 2,900 & 1970-2009 & Yes & Yes \\
20 years TW & 0.074$^{**}$ & $-$0.054 & HP & 40 | 1,450 & 1970-2009 & Yes & Yes \\
Excluding EU CP & 0.056$^{**}$ & 0.005 & HP & 40 | 2,280 & 1970-2009 & Yes & Yes \\
Excluding USSR & 0.064$^{**}$ & $-$0.006 & HP & 34 | 2,244 & 1970-2009 & Yes & Yes \\
Alternative TW & 0.081$^{***}$ & 0.014 & HP & 34 | 2,244 & 1970-1999 & Yes & Yes \\ [0.15cm]
\multicolumn{8}{l}{\textit{Alternative controls for sectoral composition}} \\ [0.1cm]
4Digits SITC & 0.058$^{**}$ & $-$0.045$^{*}$ & HP & 36 | 2,520 & 1970-2009 & Yes & Yes \\
ISIC classification & 0.059$^{**}$ & $-$0.045$^{*}$ & HP & 36 | 2,520 & 1970-2009 & Yes & Yes \\
1Digit Agg. sectors & 0.088 & $-$0.044 & HP & 38 | 1,291 & 1970-2009 & Yes & Yes \\ [0.15cm]
\multicolumn{8}{l}{\textit{Alternative indexes}} \\ [0.1cm]
$level(trade)^{ a}$ & 33,96$^{*}$ & $-$34.92 & HP & 40 | 2,900 & 1970-2009 & Yes & Yes \\
$\log(mean(trade))$ & 0.044$^{*}$ & $-$0.027 & HP & 40 | 2,900 & 1970-2009 & Yes & Yes \\
$\max \Big( \frac{T_{i\leftrightarrow j}}{GDP_i}, \frac{T_{i\leftrightarrow j}}{GDP_j} \Big) $ & 0.052$^{**}$ & $-$0.032 & HP & 40 | 2,900 & 1970-2009 & Yes & Yes \\
STAN data & 0.209$^{**}$ & $-$0.107 & HP & 20 | 760 & 1995-2014 & Yes & Yes \\
\bottomrule \\ [-0.3cm]
\multicolumn{8}{l}{\textit{Notes:} $^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01. In parenthesis: std. deviation.} \\ 
\multicolumn{8}{l}{\parbox{15cm}{$^{a}$ We provide the results using EU and USSR dummies since adding those controls substantially reduce the significance of trade in final goods.}} \\ 
\end{tabular} 
\end{table}
