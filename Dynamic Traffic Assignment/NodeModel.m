function [TurnFlows] = NodeModel(nbIncomingLinks,nbOutgoingLinks,sendingFlow,turningFractions,receivingFlow,incomingLinks_capacity)
%#codegen
%coder.inline('never');

%Node model
%
%SYNTAX
%   [flowIncomingLinks,flowOutgoingLinks] = NodeModel(incomingLinks,outgoingLinks,sendingFlow,turningFractions,receivingFlow,caps,timeInterval)
%
%DESCRIPTION
%   returns the transfer flow over a node
%
%INPUTS
%   incomingLinks: id of the incoming links
%   outgoingLinks: id of the outgoing links
%   sendingFlow: the sending flow
%   turningFractions: turning fractions at the node
%   receivingFlow: the receiving flow
%   links: list of all the links in the network
%   Each entry of the list represents one link. Each link is a structure that
%   has at least a link ID and an upstream and downstream node.
%   timeInterval: updating interval
    
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

%STATIC INITIALIZATION WITH MAX 7 IN- AND OUTGOING LINKS

%check also here dynamic initialization based on actual number of incoming
%and outgoing links

%initialize
adjustedReceivingFlow = receivingFlow; % adjusted receiving flow
competingLinks = repmat(sendingFlow>0,1,nbOutgoingLinks);
competingLinks(turningFractions==0)=0;   % competing links
activeOutLinks = any(competingLinks,1);% active outgoing links

distrFactors = turningFractions.*repmat(incomingLinks_capacity,1,nbOutgoingLinks); %oriented capacity proportional
TurnFlows = zeros(nbIncomingLinks,nbOutgoingLinks);

while any(activeOutLinks)
    alpha=adjustedReceivingFlow./sum(distrFactors.*competingLinks,1)';
    alpha(~activeOutLinks)=inf;
    [alpha_min,j]=min(alpha);
    if any(sendingFlow(competingLinks(:,j))<=alpha_min*incomingLinks_capacity(competingLinks(:,j)))
        for a=1:nbIncomingLinks
            if competingLinks(a,j)
                if sendingFlow(a)<=alpha_min*incomingLinks_capacity(a)
                    for b=1:nbOutgoingLinks
                          if activeOutLinks(b)
                            TurnFlows(a,b)=turningFractions(a,b)*sendingFlow(a);
                            adjustedReceivingFlow(b)=adjustedReceivingFlow(b)-TurnFlows(a,b);
                            competingLinks(a,b)=false;
                            if all(~competingLinks(:,b))
                                activeOutLinks(b)=false;
                            end
                        end
                    end
                end
            end
        end
    else
        for a=1:nbIncomingLinks
            if competingLinks(a,j)
                for b=1:nbOutgoingLinks
                   if activeOutLinks(b)
                        TurnFlows(a,b)=alpha_min*distrFactors(a,b);
                        adjustedReceivingFlow(b)=adjustedReceivingFlow(b)-TurnFlows(a,b);
                        competingLinks(a,b)=false;
                        if all(~competingLinks(:,b))
                            activeOutLinks(b)=false;
                        end
                    end
                end
            end
        end
        activeOutLinks(j)=false;
    end
end
end