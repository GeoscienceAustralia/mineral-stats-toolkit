# -*- coding: utf-8 -*-
"""
   Copyright 2022 Commonwealth of Australia (Geoscience Australia)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""
import os
from toolkit.cdf import CDF




wd_outputs = r'C:\tmp\outputs'

distance_array_fn = os.path.join(wd_outputs,
                                 'orogenic_gold_100ohmm_distance_array.npy')

savepath = r'C:\tmp\outputs'

# second, compute CDF
cdf_obj = CDF(distance_array_fn = distance_array_fn,
                     cdf_array_fn = 'orogenic_gold_cdf',
                     savepath=savepath,
                     binsize=2e3
                     )

cdf_obj.compute_cdf()
cdf_obj.save_cdf_array()