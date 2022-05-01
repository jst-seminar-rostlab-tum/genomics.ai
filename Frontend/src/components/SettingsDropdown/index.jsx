import React from 'react';
import { Link } from 'react-router-dom';
import Stack from '@mui/material/Stack';
import styles from './settingsDropdown.module.css';

import { useAuth } from 'shared/context/authContext';

export function SettingsDropdown() {
  const [, setUser] = useAuth();
  return (
    <Stack
      className={styles.dropdownMenu}
      spacing={3}
      direction="column"
    >
      <Link
        className={styles.dropdownLink}
        to="/sequencer/settings"
        style={{ paddingTop: '10px' }}
      >
        Profile Settings
      </Link>

      <div className={styles.divider} />

      <Link
        className={styles.dropdownLink}
        to="/"
        onClick={() => {
          setUser(null);
          localStorage.removeItem('user');
          localStorage.removeItem('jwt');
        }}
        style={{ paddingBottom: '10px' }}
      >
        Log out
      </Link>
    </Stack>
  );
}

export default SettingsDropdown;
