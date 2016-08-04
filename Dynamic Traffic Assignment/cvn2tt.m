function simTT= cvn2tt(cvn_up,cvn_down,dt,totT,links)
    
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
        [down,iun] = unique(cvn_down(l,:),'first');
        if length(down)<=1
            simTT(l,:)=links.length(l)/links.freeSpeed(l);
        else
            simTT(l,:)=max(interp1(down,timeSteps(iun),cvn_up(l,:))-dt*[0:totT],links.length(l)/links.freeSpeed(l));
            simTT(l,cvn_up(l,:)-cvn_down(l,:)<10-3)=links.length(l)/links.freeSpeed(l);
            for t=2:totT+1
                simTT(l,t)=max(simTT(l,t),simTT(l,t-1)-dt);
            end
        end
    end
end


