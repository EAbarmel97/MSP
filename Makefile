JULIA_DEPOT_PATH := $(shell pwd)/.julenv

define julia_env
	julia --project=$(JULIA_DEPOT_PATH) $(1)
endef

createprojectdir:
	@mkdir -p $(JULIA_DEPOT_PATH)

juliaenv:
	@$(call julia_env,$(ARGS))
	
instantiatenv: createprojectdir
	@cp Project.toml $(JULIA_DEPOT_PATH)/Project.toml
	@$(call julia_env,-e 'using Pkg; Pkg.resolve(); Pkg.instantiate()')

.PHONY: createprojectdir juliaenv instantiatenv