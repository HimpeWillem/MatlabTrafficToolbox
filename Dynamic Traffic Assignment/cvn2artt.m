function simTT= cvn2artt(cvn_up,cvn_down,dt,totT,links)
        
% This file is part of the matlab package for dynamic traffic assignments 
% developed by the KULeuven. 
%
% Copyright (C) 2016  Himpe Willem, Leuven, Belgium
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% More information at: http://www.mech.kuleuven.be/en/cib/traffic/downloads
% or contact: willem.himpe {@} kuleuven.be

    simTT = zeros(length(links.id),totT+1);
    %compute the simulated travel times
    timeSteps = dt*[0:1:totT];
    for l=1:length(links.length)
        [up,iun] = unique(cvn_up(l,:),'last');
        if length(up)<=1
            simTT(l,:)=links.length(l)/links.freeSpeed(l);
        else
            simTT(l,:)=max(dt*[0:totT]-interp1(up,timeSteps(iun),cvn_down(l,:)),links.length(l)/links.freeSpeed(l));
            simTT(l,cvn_up(l,:)-cvn_down(l,:)<10-3)=links.length(l)/links.freeSpeed(l);
        end
    end
end