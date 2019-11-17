/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Crypt.h
 * Author: dreamcrash
 *
 * Created on 14 de Outubro de 2016, 17:23
 */

#ifndef CRYPT_H
#define CRYPT_H

#define byte unsigned char

typedef struct IDEATest
{
    int array_rows; 
    
    byte *plain1;       // Buffer for plaintext data.
    byte *crypt1;       // Buffer for encrypted data.
    byte *plain2;       // Buffer for decrypted data.

    short *userkey;     // Key for encryption/decryption.
    int *Z;             // Encryption subkey (userkey derived).
    int *DK;            // Decryption subkey (userkey derived).
    
    
}Ideatest;

void    run             (int size, int validation, int numThreads);
int     buildTestData   (Ideatest *data);
void    calcEncryptKey  (Ideatest *data);
void    calcDecryptKey  (Ideatest *data);
int     inv             (int x);
void    Do              (Ideatest *data);
void    cipher_idea     (byte *text1, byte *text2, int *key, int text1_lenght);
void    JGFvalidate     (Ideatest *data);
#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* CRYPT_H */

