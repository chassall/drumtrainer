function bvtriggers = num2bv(condtriggers, basetrigger)
%NUM2BV Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    basetrigger = 0;
end
triggers = [];
for i = 1:length(basetrigger)
   triggers =  [triggers basetrigger(i) + condtriggers];
end

triggersFormat = 'S%3.0d';
bvtriggers = cellstr(num2str(triggers',triggersFormat))';
end

