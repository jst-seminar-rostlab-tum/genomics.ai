import React from 'react';
import styles from './header.module.css';
import ProfileImage from 'components/ProfileImage';

function Header({
  title,
  replaceRight,
}) {
  return (
    <header className={styles.header}>
      <h1 className={styles.title}>{title}</h1>
      <div className={styles.right}>
        {replaceRight == null
          ? <ProfileImage sizePixels={44} /> : replaceRight}
      </div>
    </header>
  );
}

export default Header;
