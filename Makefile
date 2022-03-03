PLATFORM  := linux
PYPI_USER := ""
PYPI_PASS := ""
DOCKER_IMAGE_ID  := ""
DOCKER_IMAGE_TAG := "latest"

setup:
	@conda env create --file environment.$(PLATFORM).yml || conda env update --file environment.$(PLATFORM).yml

build_pypi_package:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	python setup.py sdist bdist_wheel

twine_upload: build_pypi_package
	@python setup.py sdist bdist_wheel
	@twine upload \
		--repository-url https://upload.pypi.org/legacy/ \
		-u $(PYPI_USER) \
		-p $(PYPI_PASS) \
		dist/*-py3-none-any.whl

dockerhub_upload:
	@docker tag $(DOCKER_IMAGE_ID) omixlab/bambu-qsar:$(DOCKER_IMAGE_TAG)
	@docker push omixlab/bambu-qsar:$(DOCKER_IMAGE_TAG)
