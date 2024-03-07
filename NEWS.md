# scKirby 0.1.5

## New features

* Move more functions from `phenomix` --> `scKirby`:
  - `compute_cor` --> `calc_cor`
  - `compute_snn`
  - `create_DimReducObject`
  - `get_cor`
  - `has_graph`
  - `infer_graph_name` --> `infer_graph_key`
  - `normalise`
  - `yeojohnson`
  - `is_installed`
* `dict_class`
  - Add "Graph" to "matrix" class.
* New exported functions:
  - `set_obsp`
  - `to_delayedarray`

# scKirby 0.1.4

## New features

- Rename `get_data` --> `get_x` to align with `anndata` nomenclature 
and avoid conflicts with piggyback-oriented functions.
- `is_class`
  - Add new classes: "matrix_strict", "matrix_sparse", "delayed_array"
* `get_x`
    - Add new args `slot` and `assay`
* New funcs: 
  - `get_obs_names`
  - `get_var_names`
  - `mofa_to_dimreduc`
* `get_obsm`: can now handle MOFA models.
* `get_varm`: can now handle MOFA models.
* `get_obs`/ `get_var` can now take matrices.
* `is_class`:
  - Now recognises class "matrix_list"


# scKirby 0.1.3

## New features

* `anndata_to_ctd`:
  - New function to create CTD by chunking anndata obj.
* `seurat_to_anndata`
  - Add support for `anndataR` method.
  
## Bug fixes

* `map_data`
  - Change default: `as_sparse=TRUE`

# scKirby 0.1.2

## New features

* Move `scNLP::seurat_pipeline` --> `scKirby::process_seurat`
  - New subfunction `calc_specificity`

# scKirby 0.1.1

* Implemented `rworkflows` 
* Added a `NEWS.md` file to track changes to the package.
