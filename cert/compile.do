local SRC ../src
local TEST ../cert

cd `SRC'
cap do lpils
if _rc {
	cap noi do lpils
	cd `TEST'
	exit 198
}

cd `TEST'

adopath ++ `SRC'
