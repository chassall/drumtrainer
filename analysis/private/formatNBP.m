function formatNBP(h,plotColours, plotAlpha)
%FORMATNBP Formats a notBoxPlot
%   Formats colours, alpha values, line widths, and point size of the notBoxPlot.
%
%   h: handle to notBoxPlot
%   plotColours (optional): desired colours (you need to provide at least as many as
%   your conditions) as a numberOfColours X 3 matrix.
%
%   C.D. Hassall, 2018

for i = 1:length(h)
    h(i).sdPtch.FaceColor = 'none';
    h(i).sdPtch.LineStyle = 'none';
    if nargin == 1
        h(i).semPtch.FaceAlpha = 0.2;
        h(i).semPtch.FaceColor = 'none';
    elseif nargin == 2
        h(i).semPtch.FaceColor = plotColours(i,:);
        h(i).semPtch.FaceAlpha = 0.2;
    else
        h(i).semPtch.FaceColor = plotColours(i,:);
        h(i).semPtch.FaceAlpha = plotAlpha;
    end
    h(i).semPtch.LineWidth = 0.5;
    h(i).semPtch.EdgeColor = 'black';
    h(i).mu.Color = 'black';
    h(i).mu.LineWidth = 1;
    h(i).data.MarkerSize = 2;
    h(i).data.Marker = '.';
    h(i).data.MarkerFaceColor = 'none';
end

end

