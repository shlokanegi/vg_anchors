.PHONY: init sdust clean

init: sdust
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

clean:
	rm -rf bin/ third_party/sdust/sdust third_party/sdust/sdust.o
