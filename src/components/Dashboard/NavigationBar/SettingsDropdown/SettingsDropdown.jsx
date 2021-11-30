import React from 'react';
import { Link } from 'react-router-dom';
import Stack from '@mui/material/Stack';
import styles from './settingsDropdown.module.css';

export function SettingsDropdown() {
  return (
    <Stack
      className={styles.dropdownMenu}
      spacing={3}
      direction="column"
    >
      <Link
        className={styles.dropdownLink}
        to="/settings"
        style={{ paddingTop: '10px' }}
      >
        Profile Settings
      </Link>

      <Link
        className={styles.dropdownLink}
        to="/marketing"
      >
        View Profile
      </Link>

      <div className={styles.divider} />

      <Link
        className={styles.dropdownLink}
        to="/logOut"
        style={{ paddingBottom: '10px' }}
      >
        Log out
      </Link>
    </Stack>
  );
}

export default SettingsDropdown;
