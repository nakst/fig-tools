// This file is part of the Essence operating system.
// It is released under the terms of the MIT license -- see LICENSE.md.
// Written by: nakst.

struct _ArrayHeader {
	size_t length, allocated;
};

bool _ArrayEnsureAllocated(void **array, size_t minimumAllocated, size_t itemSize, uint8_t additionalHeaderBytes, EsHeap *heap);
bool _ArraySetLength(void **array, size_t newLength, size_t itemSize, uint8_t additionalHeaderBytes, EsHeap *heap);
void _ArrayDelete(void *array, uintptr_t position, size_t itemSize, size_t count);
void _ArrayDeleteSwap(void *array, uintptr_t position, size_t itemSize);
void *_ArrayInsert(void **array, const void *item, size_t itemSize, ptrdiff_t position, uint8_t additionalHeaderBytes, EsHeap *heap);
void *_ArrayInsertMany(void **array, size_t itemSize, ptrdiff_t position, size_t count, EsHeap *heap);
void _ArrayFree(void **array, size_t itemSize, EsHeap *heap);

#define ArrayHeader(array) ((_ArrayHeader *) (array) - 1)
#define ArrayLength(array) ((array) ? (ArrayHeader(array)->length) : 0)

#ifdef KERNEL
template <class T, EsHeap *heap>
#else
template <class T, EsHeap *heap = nullptr>
#endif
struct Array {
	T *array;

	T &First() { return array[0]; }
	T &Last() { return array[Length() - 1]; }
	void Delete(uintptr_t position) { _ArrayDelete(array, position, sizeof(T), 1); }
	void DeleteSwap(uintptr_t position) { _ArrayDeleteSwap(array, position, sizeof(T)); }
	void DeleteMany(uintptr_t position, size_t count) { _ArrayDelete(array, position, sizeof(T), count); }
	T *Add(T item) { return (T *) _ArrayInsert((void **) &array, &item, sizeof(T), -1, 0, heap); }
	T *Add() { return (T *) _ArrayInsert((void **) &array, nullptr, sizeof(T), -1, 0, heap); }
	T *Insert(T item, uintptr_t position) { return (T *) _ArrayInsert((void **) &array, &item, sizeof(T), position, 0, heap); }
	T *AddPointer(const T *item) { return (T *) _ArrayInsert((void **) &(array), item, sizeof(T), -1, 0, heap); }
	T *InsertPointer(const T *item, uintptr_t position) { return (T *) _ArrayInsert((void **) &array, item, sizeof(T), position, 0, heap); }
	T *InsertMany(uintptr_t position, size_t count) { return (T *) _ArrayInsertMany((void **) &array, sizeof(T), position, count, heap); }
	bool SetLength(size_t length) { return _ArraySetLength((void **) &array, length, sizeof(T), 0, heap); }
	void Free() { _ArrayFree((void **) &array, sizeof(T), heap); }
	T Pop() { T t = Last(); Delete(Length() - 1); return t; }

	__attribute__((no_instrument_function)) 
	size_t Length() { 
		return array ? ArrayHeader(array)->length : 0; 
	}

	__attribute__((no_instrument_function)) 
	T &operator[](uintptr_t index) { 
#ifdef DEBUG_BUILD
		EsAssert(index < Length());
#endif
		return array[index]; 
	}

	intptr_t Find(T item, bool failIfNotFound) {
		for (uintptr_t i = 0; i < Length(); i++) {
			if (array[i] == item) {
				return i;
			}
		}

		if (failIfNotFound) EsPanic("Array::Find - Item not found in %x.\n", this);
		return -1;
	}

	bool FindAndDelete(T item, bool failIfNotFound) {
		intptr_t index = Find(item, failIfNotFound);
		if (index == -1) return false;
		Delete(index);
		return true;
	}

	bool FindAndDeleteSwap(T item, bool failIfNotFound) { 
		intptr_t index = Find(item, failIfNotFound);
		if (index == -1) return false;
		DeleteSwap(index);
		return true;
	}

	void AddFast(T item) { 
		if (!array) { Add(item); return; }
		_ArrayHeader *header = ArrayHeader(array);
		if (header->length == header->allocated) { Add(item); return; }
		array[header->length++] = item;
	}
};

#ifndef ARRAY_DEFINITIONS_ONLY

bool _ArrayMaybeInitialise(void **array, size_t itemSize, EsHeap *heap) {
	if (*array) return true;
	size_t newLength = 4;
	_ArrayHeader *header = (_ArrayHeader *) EsHeapAllocate(sizeof(_ArrayHeader) + itemSize * newLength, true, heap);
	if (!header) return false;
	header->length = 0;
	header->allocated = newLength;
	*array = header + 1;
	return true;
}

bool _ArrayEnsureAllocated(void **array, size_t minimumAllocated, size_t itemSize, uint8_t additionalHeaderBytes, EsHeap *heap) {
	if (!_ArrayMaybeInitialise(array, itemSize, heap)) {
		return false;
	}

	_ArrayHeader *oldHeader = ArrayHeader(*array);

	if (oldHeader->allocated >= minimumAllocated) {
		return true;
	}

	_ArrayHeader *newHeader = (_ArrayHeader *) EsHeapReallocate((uint8_t *) oldHeader - additionalHeaderBytes, 
			sizeof(_ArrayHeader) + additionalHeaderBytes + itemSize * minimumAllocated, false, heap);

	if (!newHeader) {
		return false;
	}

	newHeader->allocated = minimumAllocated;
	*array = (uint8_t *) (newHeader + 1) + additionalHeaderBytes;
	return true;
}

bool _ArraySetLength(void **array, size_t newLength, size_t itemSize, uint8_t additionalHeaderBytes, EsHeap *heap) {
	if (!_ArrayMaybeInitialise(array, itemSize, heap)) {
		return false;
	}

	_ArrayHeader *header = ArrayHeader(*array);

	if (header->allocated >= newLength) {
		header->length = newLength;
		return true;
	}

	if (!_ArrayEnsureAllocated(array, header->allocated * 2 > newLength ? header->allocated * 2 : newLength + 16, itemSize, additionalHeaderBytes, heap)) {
		return false;
	}

	header = ArrayHeader(*array);
	header->length = newLength;
	return true;
}

void _ArrayDelete(void *array, uintptr_t position, size_t itemSize, size_t count) {
	if (!count) return;
	size_t oldArrayLength = ArrayLength(array);
	if (position >= oldArrayLength) EsPanic("_ArrayDelete - Position out of bounds (%d/%d).\n", position, oldArrayLength);
	if (count > oldArrayLength - position) EsPanic("_ArrayDelete - Count out of bounds (%d/%d/%d).\n", position, count, oldArrayLength);
	ArrayHeader(array)->length = oldArrayLength - count;
	uint8_t *data = (uint8_t *) array;
	EsMemoryMove(data + itemSize * (position + count), data + itemSize * oldArrayLength, ES_MEMORY_MOVE_BACKWARDS itemSize * count, false);
}

void _ArrayDeleteSwap(void *array, uintptr_t position, size_t itemSize) {
	size_t oldArrayLength = ArrayLength(array);
	if (position >= oldArrayLength) EsPanic("_ArrayDeleteSwap - Position out of bounds (%d/%d).\n", position, oldArrayLength);
	ArrayHeader(array)->length = oldArrayLength - 1;
	uint8_t *data = (uint8_t *) array;
	EsMemoryCopy(data + itemSize * position, data + itemSize * ArrayLength(array), itemSize);
}

void *_ArrayInsert(void **array, const void *item, size_t itemSize, ptrdiff_t position, uint8_t additionalHeaderBytes, EsHeap *heap) {
	size_t oldArrayLength = ArrayLength(*array);
	if (position == -1) position = oldArrayLength;
	if (position < 0 || (size_t) position > oldArrayLength) EsPanic("_ArrayInsert - Position out of bounds (%d/%d).\n", position, oldArrayLength);
	if (!_ArraySetLength(array, oldArrayLength + 1, itemSize, additionalHeaderBytes, heap)) return nullptr;
	uint8_t *data = (uint8_t *) *array;
	EsMemoryMove(data + itemSize * position, data + itemSize * oldArrayLength, itemSize, false);
	if (item) EsMemoryCopy(data + itemSize * position, item, itemSize);
	else EsMemoryZero(data + itemSize * position, itemSize);
	return data + itemSize * position;
}

void *_ArrayInsertMany(void **array, size_t itemSize, ptrdiff_t position, size_t insertCount, EsHeap *heap) {
	size_t oldArrayLength = ArrayLength(*array);
	if (position == -1) position = oldArrayLength;
	if (position < 0 || (size_t) position > oldArrayLength) EsPanic("_ArrayInsertMany - Position out of bounds (%d/%d).\n", position, oldArrayLength);
	if (!_ArraySetLength(array, oldArrayLength + insertCount, itemSize, 0, heap)) return nullptr;
	uint8_t *data = (uint8_t *) *array;
	EsMemoryMove(data + itemSize * position, data + itemSize * oldArrayLength, itemSize * insertCount, false);
	return data + itemSize * position;
}

void _ArrayFree(void **array, size_t itemSize, EsHeap *heap) {
	if (!(*array)) return;
	EsHeapFree(ArrayHeader(*array), sizeof(_ArrayHeader) + itemSize * ArrayHeader(*array)->allocated, heap);
	*array = nullptr;
}

#endif
