function exportDemoFigure(fig, basePath, preservePageSize)
%EXPORTDEMOFIGURE Export vector PDF and a 300-dpi PNG preview.

if nargin < 3
    preservePageSize = false;
end

outDir = fileparts(basePath);
if ~isfolder(outDir)
    mkdir(outDir);
end

set(fig, 'Color', 'white', 'Renderer', 'painters', ...
    'InvertHardcopy', 'off');
if preservePageSize
    set(fig, 'PaperPositionMode', 'manual');
    print(fig, [basePath '.pdf'], '-dpdf', '-painters');
else
    set(fig, 'PaperPositionMode', 'auto');
    exportgraphics(fig, [basePath '.pdf'], 'ContentType', 'vector');
end
exportgraphics(fig, [basePath '.png'], 'Resolution', 300);
end
