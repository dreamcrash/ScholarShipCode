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
    int p_array_rows;
    int ref_p_array_rows;
    int rem_p_array_rows;
    
    // Master global data
    byte *plain1;       // Buffer for plaintext data.
    byte *crypt1;       // Buffer for encrypted data.
    byte *plain2;       // Buffer for decrypted data.
    
    // Slave/ Master copies
    byte *p_plain1;       // Buffer for plaintext data.
    byte *p_crypt1;       // Buffer for encrypted data.
    byte *p_plain2;       // Buffer for decrypted data.

    short *userkey;     // Key for encryption/decryption.
    int *Z;             // Encryption subkey (userkey derived).
    int *DK;            // Decryption subkey (userkey derived).
    
    
}Ideatest;

void    run             (int size, int validation);
int     buildTestData   (Ideatest *data, int rank, int total_process);
void    calcEncryptKey  (Ideatest *data);
void    calcDecryptKey  (Ideatest *data);
int     inv             (int x);
void    Do              (Ideatest *data, int total_process);
void    cipher_idea     (byte *text1, byte *text2, int *key, int text1_lenght);
void    freeIdeatest    (Ideatest **data);
void    JGFvalidate     (Ideatest *data);
#ifdef __cplusplus
extern "C" {
#endif




#ifdef __cplusplus
}
#endif

#endif /* CRYPT_H */

