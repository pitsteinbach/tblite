#pragma once
#include "tblite/macros.h"

// Post Processing Container
typedef struct _tblite_post_processing* tblite_post_processing;

/// Construct post processing container
///
/// @param charptr: String of the post processing desired
/// @return New post processing instance
TBLITE_API_ENTRY tblite_post_processing TBLITE_API_CALL
tblite_new_post_processing(char* charptr);

/// Delete calculator
///
/// @param post_proc: Post Processing instance
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_post_processing(tblite_post_processing* post_proc);