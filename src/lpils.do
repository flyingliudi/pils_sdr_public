mata: mata clear

do pils.mata
do st_pils_obj.mata

mata:
mata mlib create lpils, dir(.) replace
mata mlib add lpils *()
mata mlib index
end
