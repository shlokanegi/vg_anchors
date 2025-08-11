.PHONY: init sdust clean bdsg

init: sdust bdsg
	pip install -e .

sdust:
	@if [ ! -f bin/sdust ]; then \
		echo "--- Compiling sdust ---"; \
		git submodule update --init --recursive; \
		$(MAKE) -C third_party/sdust; \
		mkdir -p bin; \
		mv third_party/sdust/sdust bin/; \
	else \
		echo "--- sdust is already compiled ---"; \
	fi

bdsg:
	@echo "--- Installing bdsg dependency ---"; \
	cd libbdsg && \
	rm -rf build/ dist/ *.egg-info/ && \
	pip install -e . && \
	cd ..

clean:
	rm -rf bin/ third_party/sdust/sdust third_party/sdust/sdust.o
	cd libbdsg && rm -rf build/ dist/ *.egg-info/ && cd ..
