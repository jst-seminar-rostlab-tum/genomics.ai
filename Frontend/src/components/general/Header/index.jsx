import React from 'react';
import styles from './header.module.css';
import { Link } from 'react-router-dom';
import ProfileImage from 'components/ProfileImage';
import Stack from '@mui/material/Stack';

function Header({
  title,
  rightOfTitle,
  replaceRight,
  loggedIn = true,
}) {
  return (
    <header className={styles.header}>
      <div className={styles.titleRow}>
        <h1 className={styles.title}>{title}</h1>
        {rightOfTitle}
      </div>
      {/* render avatar and settings page if logged in */}
      {loggedIn && (
      <Stack direction="row" spacing={4} alignItems="flex-end">
        {replaceRight}
        <Link to="/sequencer/settings">
          <ProfileImage sizePixels={42} />
        </Link>
      </Stack>
      )}
    </header>
  );
}

export default Header;
