# ensemble-scoring
### Analysis tools for PISM parameter ensemble with full-factorial sampling

These tools enable for the parameter ensemble analysis as described in

- Albrecht et al., "Glacial cycles simulation of the Antarctic Ice Sheet with PISM – Part 2: Parameter ensemble analysis", The Cryosphere, (2019)

and can be applied to PISM results in

- Albrecht, Torsten (2019): PISM parameter ensemble analysis of Antarctic Ice Sheet glacial cycle simulations, PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.909728


The python-based scoring scheme with respect to modern and paleo data is performed with 'python run.py'. It is based on [Pollard et al., 2016](https://doi.org/10.5194/gmd-9-1697-2016) and [Briggs et al., 2014](https://doi.org/10.1016/j.quascirev.2014.09.003).. The ensemble analysis calculates misfits to the paleo constraint database [AntICEdat (Briggs & Tarasov, 2013)](http://dx.doi.org/10.1016/j.quascirev.2012.11.021) and to [RAISED Consortium (2014)](https://doi.org/10.1016/j.quascirev.2014.06.025) as well as to modern ice geometry from [Bedmap2 (Fretwell et al., 2013)](https://doi.org/10.5194/tc-7-375-2013), [ice speed (Rignot et al., 2011)](https://doi.org/10.1126/science.1208336) an [GPS (Whitehouse et al., 2011)](http://dx.doi.org/10.1111/j.1365-246X.2012.05557.x). 


#### Key references:

- Briggs, R. D., & Tarasov, L. (2013). How to evaluate model-derived deglaciation chronologies: a case study using Antarctica. Quaternary Science Reviews, 63, 109-127.

- Briggs, R. D., Pollard, D., & Tarasov, L. (2014). A data-constrained large ensemble analysis of Antarctic evolution since the Eemian. Quaternary Science Reviews, 103, 91-115.

- Pollard, D., Chang, W., Haran, M., Applegate, P., and DeConto, R. (2016): Large ensemble modeling of the last deglacial retreat of the West Antarctic Ice Sheet: comparison of simple and advanced statistical techniques, Geosci. Model Dev., 9, 1697–1723


#### License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program (LICENSE). If not, see https://www.gnu.org/licenses/.



#### PISM:
PISM is the open-source Parallel Ice Sheet Model developed mainly at UAF, USA and PIK, Germany. See documentation in http://www.pism-docs.org. PISM output is in netCDF format.



#### Technical information
You may need to 'pip install scikit-fmm' or create your virtual conda environment.

[Pickle](https://wiki.python.org/moin/UsingPickle) is used to store aggregated data.

Edit config.py to adjust pathes and settings.



#### Keywords: 
Ensemble Analysis, Antarctic Ice Sheet, glacial cycles, PISM, Parallel Ice Sheet Model, Glacial Isostatic Adjustment



