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
	@echo "--- Compiling and installing libbdsg dependency ---"; \
	cd libbdsg && \
	mkdir -p build && \
	cd build && \
	cmake .. && \
	$(MAKE) -j8 && \
	cd .. && \
	pip install .

clean:
	rm -rf bin/ third_party/sdust/sdust third_party/sdust/sdust.o
	rm -rf libbdsg/build libbdsg/dist libbdsg/*.egg-info
