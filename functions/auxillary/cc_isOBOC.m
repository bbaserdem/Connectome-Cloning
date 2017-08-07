function [ map ] = cc_isOBOC( marker )
% CC_ISOBOC Returns cell-ID map if barcodes are OBOC, and false otherwise.
%   mapping is such that answer(cell_number) will give the cell as
%   identified by barcode.
%   Requires the marker matrix
% This makes it so that both cases where;
%   -Parent connectome is smaller than the tabula rasa one
% 	-Parent connectome is the same size as the tabula rasa one, but has 
%   unconnected cells
%   -Parent connectome is larger than tabula-rasa, but has empty cells to 
%   accomodate an OBOC solution on tabula-rasa network
%       -> If the parent is larger, but connected is as much as child, no
%       emptiness will occur
%       -> If the parent is same size, but connected is less than child;
%       some barcodes will be empty and they will show up in the difference
%       betwee {1,2,...,child}.
%       -> If the parent is larger but connected is less than child;
%           -> Worst case is all connected ID's are less than child_size.
%           Then unconnected will appear within the range [1, child]
%           sufficiently that random assigning will handle it
%           -> If disconnected barcodes have larger ID# than child_size,
%           they will increase the number of available free assigned ID's.
%           The difference [1,...,child_size]\[id'd_barcodes] will be
%           larger than empty cells in child
% Marker matrix of form (barcode ID, cell#)

% Initialize
map = false;

% Colsum should be equal to colmax
if any( max(marker) ~= sum(marker) )
    % Ensures that all barcodes in a cell is a single barcode type
    return
% Rowsum should be equal to rowmax
elseif any( max(marker,[],2) ~= sum(marker,2) )
    % Ensures that a barcode type is only present in a single cell
    return
end

% Not interested in how many barcodes now, but which barcode ID is present
[~,map] = max( marker );
% Remove the false 1 ID's of empty cells by O'ing their identification ID
map( sum(marker) == 0 ) = 0;
% Check which ones are missing ID's
miss = find( map==0 );
% Usable ID numbers for identification 
uuid = setdiff( 1:size(marker,2), map(map~=0) );
% Identify empty cells in order
map( miss ) = uuid( 1:length(miss) );

end
