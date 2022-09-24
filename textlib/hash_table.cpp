// This file is part of the Essence operating system.
// It is released under the terms of the MIT license -- see LICENSE.md.
// Written by: nakst.

// TODO Test deleting values in hash stores with variable length keys.

//////////////////////////////////////////

struct HashTableKey {
	union {
		uintptr_t shortKey;
		void *longKey;
	};

	uint32_t longKeyBytes;
	uint16_t longKeyHash16;
	uint8_t used, tombstone;
};

struct HashTableSlot {
	HashTableKey key;
	void *value;
};

struct HashTable {
	size_t itemCount;
	size_t slotCount;
	HashTableSlot *slots;
	uint8_t *storage;
};

struct HashStoreOptions {
	size_t keyBytes;
	size_t valueBytes;
};

uint64_t HashTableFNV1a(const void *key, size_t keyBytes) {
	uint64_t hash = 0xCBF29CE484222325;

	for (uintptr_t i = 0; i < keyBytes; i++) {
		hash = (hash ^ ((uint8_t *) key)[i]) * 0x100000001B3;
	}

	return hash;
}

HashTableSlot *_HashTableGetSlot(HashTable *table, HashTableKey *key, bool useLongKeys) {
	uint64_t hash;

	if (useLongKeys) {
		hash = HashTableFNV1a(key->longKey, key->longKeyBytes);
		key->longKeyHash16 = (uint16_t) (hash >> 48);
	} else {
		hash = HashTableFNV1a(&key->shortKey, sizeof(key->shortKey));
	}

	uintptr_t slot = hash & (table->slotCount - 1);
	HashTableSlot *firstTombstone = nullptr;

	for (uintptr_t probe = 0; probe < table->slotCount; probe++) {
		uintptr_t index = slot + probe;
		if (index >= table->slotCount) index -= table->slotCount;
		HashTableSlot *s = table->slots + index;

		if (!s->key.used && !s->key.tombstone) {
			return firstTombstone ?: s;
		}

		if (s->key.tombstone) {
			if (!firstTombstone) {
				firstTombstone = s;
			}

			continue;
		}

		if (useLongKeys) {
			if (s->key.longKeyHash16 == key->longKeyHash16 
					&& s->key.longKeyBytes == key->longKeyBytes 
					&& 0 == EsMemoryCompare(s->key.longKey, key->longKey, key->longKeyBytes)) {
				return s;
			}
		} else {
			if (s->key.shortKey == key->shortKey) {
				return s;
			}
		}
	}

	EsAssert(firstTombstone);
	return firstTombstone;
}

HashTableSlot *HashTableGetSlot(HashTable *table, HashTableKey key, bool useLongKeys) {
	if (!table->itemCount) return nullptr;
	HashTableSlot *slot = _HashTableGetSlot(table, &key, useLongKeys);
	return slot->key.used ? slot : nullptr;
}

void *HashTableGet(HashTable *table, HashTableKey key, bool useLongKeys) {
	HashTableSlot *slot = HashTableGetSlot(table, key, useLongKeys);
	return slot ? slot->value : nullptr;
}

bool _HashTableResize(HashTable *table, size_t newSlotCount, bool useLongKeys);

bool _HashTablePutSlot(HashTable *table, HashTableSlot *slot, HashTableKey key, void *value, bool useLongKeys, bool duplicateLongKeys) {
	if (slot->key.used) {
		if (!duplicateLongKeys) {
			slot->key = key;
			slot->key.tombstone = false;
			slot->key.used = true;
		}

		slot->value = value;
	} else {
		if (useLongKeys && duplicateLongKeys) {
			void *keyCopy = EsHeapAllocate(key.longKeyBytes, false);
			if (!keyCopy) return false;
			EsMemoryCopy(keyCopy, key.longKey, key.longKeyBytes);
			key.longKey = keyCopy;
		}

		slot->key = key;
		slot->key.tombstone = false;
		slot->key.used = true;
		slot->value = value;
		table->itemCount++;
	}

	return true;
}

bool _HashTableEnsureSpaceAvailable(HashTable *table, bool useLongKeys) {
	if (!table->slotCount) {
		table->slotCount = 16;
		table->slots = (HashTableSlot *) EsHeapAllocate(sizeof(HashTableSlot) * table->slotCount, true);

		if (!table->slots) {
			return false;
		}
	}

	if (table->itemCount > table->slotCount / 2) {
		if (!_HashTableResize(table, table->slotCount * 2, useLongKeys)) {
			return false;
		}
	}

	return true;
}

bool HashTablePut(HashTable *table, HashTableKey key, void *value, bool useLongKeys, bool duplicateLongKeys = true) {
	if (!_HashTableEnsureSpaceAvailable(table, useLongKeys)) return false;
	HashTableSlot *slot = _HashTableGetSlot(table, &key, useLongKeys);
	return _HashTablePutSlot(table, slot, key, value, useLongKeys, duplicateLongKeys);
}

bool _HashTableResize(HashTable *table, size_t newSlotCount, bool useLongKeys) {
	HashTable copy = *table;
	copy.slotCount = newSlotCount;
	copy.slots = (HashTableSlot *) EsHeapAllocate(newSlotCount * sizeof(HashTableSlot), true);
	copy.itemCount = 0;

	if (!copy.slots) {
		return false;
	}

	for (uintptr_t i = 0; i < table->slotCount; i++) {
		if (table->slots[i].key.used) {
			bool success = HashTablePut(&copy, table->slots[i].key, table->slots[i].value, useLongKeys, false);
			EsAssert(success);
		}
	}

	EsAssert(table->itemCount == copy.itemCount);
	EsHeapFree(table->slots);
	*table = copy;
	return true;
}

bool HashTableDelete(HashTable *table, HashTableKey key, bool useLongKeys) {
	HashTableSlot *slot = HashTableGetSlot(table, key, useLongKeys);
	if (!slot) return false;
	slot->key.used = false;
	slot->key.tombstone = true;
	table->itemCount--;

	if (table->itemCount < table->slotCount / 8 && table->itemCount > 16) {
		_HashTableResize(table, table->slotCount / 2, useLongKeys);
	}

	return true;
}

void HashTableFree(HashTable *table, bool freeLongKeys) {
	if (freeLongKeys) {
		for (uintptr_t i = 0; i < table->slotCount; i++) {
			if (table->slots[i].key.used) {
				table->itemCount--;
				EsHeapFree(table->slots[i].key.longKey);
			}
		}

		EsAssert(!table->itemCount);
	}

	EsHeapFree(table->slots);
	if (table->storage) EsHeapFree(table->storage - sizeof(size_t));
	table->itemCount = 0;
	table->slotCount = 0;
	table->slots = nullptr;
	table->storage = nullptr;
}

//////////////////////////////////////////

void *HashTableGetLong(HashTable *table, const void *key, size_t keyBytes) {
	HashTableKey k;
	k.longKey = (void *) key, k.longKeyBytes = keyBytes;
	return HashTableGet(table, k, true);
}

void *HashTableGetShort(HashTable *table, uintptr_t key) {
	HashTableKey k;
	k.shortKey = key;
	return HashTableGet(table, k, false);
}

bool HashTablePutLong(HashTable *table, const void *key, size_t keyBytes, void *value, bool duplicateKey = true) {
	HashTableKey k;
	k.longKey = (void *) key, k.longKeyBytes = keyBytes;
	return HashTablePut(table, k, value, true, duplicateKey);
}

bool HashTablePutShort(HashTable *table, uintptr_t key, void *value) {
	HashTableKey k;
	k.shortKey = key;
	return HashTablePut(table, k, value, false);
}

bool HashTableDeleteLong(HashTable *table, const void *key, size_t keyBytes) {
	HashTableKey k;
	k.longKey = (void *) key, k.longKeyBytes = keyBytes;
	return HashTableDelete(table, k, true);
}

bool HashTableDeleteShort(HashTable *table, uintptr_t key) {
	HashTableKey k;
	k.shortKey = key;
	return HashTableDelete(table, k, false);
}

//////////////////////////////////////////

void *HashStoreGet(HashTable *table, HashStoreOptions *options, const void *key, size_t keyBytes = 0) {
	return HashTableGetLong(table, key, options->keyBytes ?: keyBytes);
}

void *HashStorePut(HashTable *table, HashStoreOptions *options, const void *key, size_t keyBytes = 0) {
	if (!_HashTableEnsureSpaceAvailable(table, true)) return nullptr;
	bool variableLengthKeys = keyBytes;
	if (!keyBytes) keyBytes = options->keyBytes;

	HashTableKey k;
	k.longKey = (void *) key, k.longKeyBytes = keyBytes;
	HashTableSlot *slot = _HashTableGetSlot(table, &k, true);

	if (!slot->key.used) {
		if (!table->storage) {
			table->storage = (uint8_t *) EsHeapAllocate(
					(options->keyBytes + options->valueBytes) * 16 + sizeof(size_t), true);
			if (!table->storage) return nullptr;
			table->storage += sizeof(size_t);
			((size_t *) table->storage)[-1] = 16;
		}

		size_t storageAllocated = ((size_t *) table->storage)[-1];

		if (table->itemCount == storageAllocated) {
			size_t oldSize = (options->keyBytes + options->valueBytes) * storageAllocated;
			size_t newSize = oldSize * 2;
			uint8_t *newStorage = (uint8_t *) EsHeapReallocate(table->storage - sizeof(size_t), newSize + sizeof(size_t), false, nullptr);
			if (!newStorage) return nullptr;
			newStorage += sizeof(size_t);
			EsMemoryZero(newStorage + oldSize, newSize - oldSize);

			if (newStorage != table->storage) {
				for (uintptr_t i = 0; i < table->slotCount; i++) {
					if (!table->slots[i].key.used) continue;

					if (!variableLengthKeys) {
						table->slots[i].key.longKey = (uint8_t *) table->slots[i].key.longKey - table->storage + newStorage;
					}

					table->slots[i].value = (uint8_t *) table->slots[i].value - table->storage + newStorage;
				}
			}

			table->storage = newStorage;
			((size_t *) table->storage)[-1] = storageAllocated * 2;
		}

		uintptr_t index = table->itemCount;
		uint8_t *value = table->storage + index * (options->keyBytes + options->valueBytes);
		uint8_t *keyDestination = value + options->valueBytes;

		if (!variableLengthKeys) {
			EsMemoryCopy(keyDestination, key, keyBytes);
			k.longKey = keyDestination;
		}

		if (!_HashTablePutSlot(table, slot, k, value, true, variableLengthKeys)) {
			return nullptr;
		}

		if (variableLengthKeys) {
			EsAssert(options->keyBytes == sizeof(void *));
			EsMemoryCopy(keyDestination, slot->key.longKey, sizeof(void *));
		}
	}

	return slot->value;
}

bool HashStoreDelete(HashTable *table, HashStoreOptions *options, const void *key, size_t keyBytes = 0) {
	if (!keyBytes) keyBytes = options->keyBytes;
	uint8_t *value = (uint8_t *) HashTableGetLong(table, key, keyBytes);
	if (!value) return false;
	bool success = HashTableDeleteLong(table, key, keyBytes);
	EsAssert(success);
	uint8_t *end = table->storage + (options->keyBytes + options->valueBytes) * table->itemCount;
	if (value == end) return true;
	EsMemoryCopy(value, end, options->keyBytes + options->valueBytes);
	HashTableKey k;
	k.longKey = value + options->valueBytes, k.longKeyBytes = keyBytes;
	success = HashTablePut(table, k, value, true, false);
	EsAssert(success);
	return true;
}

//////////////////////////////////////////

template <class K /* set to char for long keys */, class V>
struct HashStore {
	HashTable table;

	inline size_t Count() {
		return table.itemCount;
	}

	inline V &operator[](uintptr_t index) {
		uintptr_t keyBytes = sizeof(K) == 1 ? sizeof(void *) : sizeof(K);
		return *(V *) (table.storage + (sizeof(V) + keyBytes) * index);
	}

	inline K KeyAtIndex(uintptr_t index) {
		uintptr_t keyBytes = sizeof(K) == 1 ? sizeof(void *) : sizeof(K);
		return *(K *) (table.storage + (sizeof(V) + keyBytes) * index + sizeof(V));
	}

	inline V *Get(const K *k, size_t bytes = 0) {
		HashStoreOptions options = { sizeof(K), sizeof(V) };
		if (options.keyBytes == 1) options.keyBytes = 0;
		return (V *) HashStoreGet(&table, &options, k, bytes);
	}

	inline V Get1(const K *k, size_t bytes = 0) {
		HashStoreOptions options = { sizeof(K), sizeof(V) };
		if (options.keyBytes == 1) options.keyBytes = 0;
		V *v = (V *) HashStoreGet(&table, &options, k, bytes);
		if (v) return *v;
		else return {};
	}

	inline V *Put(const K *k, size_t bytes = 0) {
		HashStoreOptions options = { sizeof(K), sizeof(V) };
		if (options.keyBytes == 1) options.keyBytes = sizeof(void *);
		return (V *) HashStorePut(&table, &options, k, bytes);
	}

	inline bool Delete(const K *k, size_t bytes = 0) {
		HashStoreOptions options = { sizeof(K), sizeof(V) };
		if (options.keyBytes == 1) options.keyBytes = sizeof(void *);
		return HashStoreDelete(&table, &options, k, bytes);
	}

	inline void Free() {
		HashTableFree(&table, sizeof(K) == 1);
	}
};

//////////////////////////////////////////

#if 0

#include <stdio.h>
#include "stb_ds.h"

struct TestItem {
	int key, value;
};

TestItem *testItems;
HashTable table;

void TestGet(int key) {
	int value1 = hmget(testItems, key);
	int value2 = (intptr_t) HashTableGetShort(&table, key);
	assert(value1 == value2);
}

void TestCompare() {
	for (int i = 0; i < 256; i++) {
		TestGet(i);
	}

	assert(hmlenu(testItems) == table.itemCount);
}

void TestPut(int key, int value) {
	hmput(testItems, key, value);
	HashTablePutShort(&table, key, (void *) (intptr_t) value);
	TestCompare();
}

void TestDelete(int key) {
	TestGet(key);
	hmdel(testItems, key);
	HashTableDeleteShort(&table, key);
	TestCompare();
}

int main() {
	{
		HashTable store = {};
		HashStoreOptions options = { .keyBytes = sizeof(uint64_t), .valueBytes = sizeof(float) };
		uint64_t key;

		key = 1;
		*(float *) HashStorePut(&store, &options, &key) = 5.3f;
		assert(*(float *) HashStoreGet(&store, &options, &key) == 5.3f);
		key = 2;
		assert(!HashStoreGet(&store, &options, &key));

		for (uintptr_t i = 0; i < 100; i++) {
			key = i;
			*(float *) HashStorePut(&store, &options, &key) = i;
		}

		for (uintptr_t i = 0; i < 100; i++) {
			key = i;
			assert(*(float *) HashStoreGet(&store, &options, &key) == i);
		}

		for (uintptr_t i = 1; i < 100; i += 2) {
			key = i;
			assert(HashStoreDelete(&store, &options, &key));
		}

		for (uintptr_t i = 0; i < 100; i += 2) {
			key = i;
			assert(*(float *) HashStoreGet(&store, &options, &key) == i);
		}

		for (uintptr_t i = 1; i < 100; i += 2) {
			key = i;
			assert(!HashStoreGet(&store, &options, &key));
		}

		for (uintptr_t i = 0; i < 100; i += 2) {
			key = i;
			assert(HashStoreDelete(&store, &options, &key));
		}

		for (uintptr_t i = 0; i < 100; i++) {
			key = i;
			assert(!HashStoreGet(&store, &options, &key));
		}

		HashTableFree(&store, true);
	}

#if 0
	for (int i = 0; i < 127; i++) {
		TestPut(i, rand() % 255);
	}

	for (int i = 0; i < 127; i++) {
		TestDelete(i);
	}

	for (uintptr_t i = 0; i < 1000; i++) {
		if (rand() & 1) {
			for (uintptr_t i = 0; i < 100; i++) {
				TestPut(rand() % 255, rand() % 255);
			}
		} else {
			for (uintptr_t i = 0; i < 100; i++) {
				TestDelete(rand() % 127);
			}
		}
	}

	hmfree(testItems);
	HashTableFree(&table, false);
#endif

	return 0;
}

#endif
