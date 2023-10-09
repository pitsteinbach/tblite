#pragma once
#include "tblite/macros.h"

// Double Dictionary Container
typedef struct _tblite_double_dictionary* tblite_double_dictionary;

/// Get number of entries
///
/// @param dict: Double dictionary instance
TBLITE_API_ENTRY int TBLITE_API_CALL
tblite_get_n_entries_dict(tblite_double_dictionary dict);

/// Get the array associated with an entry by index, together with its dimensions
///
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param array: Array associated to the entry addressed by the index
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_entry_index(tblite_double_dictionary dict,
                            const int* index,
                            double* array);

/// Get the array associated with an entry by index, together with its dimensions
///
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param dim1: 1st dimension of the associated tensor to the index
/// @param dim2: 2nd dimension of the associated tensor to the index
/// @param dim3: 3rd dimension of the associated tensor to the index
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_array_size_index(tblite_double_dictionary dict,
                            const int* index,
                            int* dim1,
                            int* dim2,
                            int* dim3);

/// Get label of an entry by index
///
/// @param dict: Double dictionary instance
/// @param index: Index of the entry for which to retrieve the label
/// @param label: Label which is retrieved
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_get_label_entry_index(tblite_double_dictionary dict,
                            const int* index,
                            char* label);

/// Delete dictionary
///
/// @param dict: Double dictionary instance
TBLITE_API_ENTRY void TBLITE_API_CALL
tblite_delete_double_dictionary(tblite_double_dictionary* dict);