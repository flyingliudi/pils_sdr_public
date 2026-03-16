*! version 1.0.0  06feb2026
program pils_obj
	version 16.0

	gettoken subcmd tmp : 0

	if (`"`tmp'"' != "") {
		di as err "subcommand {bf:`tmp'} not allowed"
		exit 198
	}

	if (`"`subcmd'"' == "create") {
		CreateObject
	}
	else if (`"`subcmd'"' == "rm") {
		RmObject
	}
	else if (`"`subcmd'"' == "getname") {
		GetName	
		c_local OBJ `s(pils_objname)'
	}
end

					//----------------------------//
					// Create object
					//----------------------------//
program CreateObject
	mata: st_pils_create_object()
end
					//----------------------------//
					// rm object
					//----------------------------//
program RmObject
	mata: st_pils_clean_object()
end
					//----------------------------//
					// get name
					//----------------------------//
program GetName
	mata: st_pils_get_object()
end
