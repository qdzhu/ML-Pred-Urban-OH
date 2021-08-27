function [ tuv_eqn ] = wrf_tuv_eqn_mapper( wrf_j )
%[ TUV_EQN ] WRF_TUV_EQN_MAPPER( WRF_J )
%   Returns the equation from TUV that corresponds to the desired WRF-Chem
%   photolysis reaction. WRF_J must be the name of the WRF-Chem KPP
%   photolysis constant, that is, if in the .eqn file, it is listed as
%   'j(Pj_no2)' pass in j(Pj_no2) or Pj_no2

[s,e] = regexp(wrf_j, '\(.*\)');
if ~isempty(s)
    wrf_j = wrf_j(s+1:e-1);
end

switch wrf_j
    case 'Pj_no2'
        tuv_eqn = 'NO2 -> NO + O(3P)';
end

end

