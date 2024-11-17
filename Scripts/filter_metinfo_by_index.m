function met_info = filter_metinfo_by_index(met_info, filter_index)

met_info.MZ = met_info.MZ(filter_index,:);
met_info.RT = met_info.RT(filter_index,:);
met_info.CompoundID = met_info.CompoundID(filter_index,:);
met_info.CompoundName = met_info.CompoundName(filter_index,:);
met_info.MetaboliteFilter = met_info.MetaboliteFilter(filter_index,:);
met_info.SumGITclusters = met_info.SumGITclusters(filter_index,:);
