function handle = subplot_handle( row, col )
	%% Initialisation of subplot function that only requires row and col in the following
	% Usage: - h = SUBPLOT_HANDLE( 2,3 ) generates subplot handle for 2 rows and 3 cols
	%		 - h(1,2) then will generate a subplot in row 1 and col 2 of this figure (even if you use a
	%			different figure in the meanwhile)
	%		 - h(1) will create a subplot with original indexing in row 1 col 1, and
	%			iterates thence through all further subplots
	% 
	% See also SUBPLOT

	cf = figure();
	handle = @create_handle;
	function create_handle( rowp, colp )
		figure( cf )
		if nargin == 1
			subplot( row, col, rowp );
		elseif nargin == 2
			subplot( row, col, (rowp-1)*col+colp );
		end
	end
	
end