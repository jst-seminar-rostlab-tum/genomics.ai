import React, { useState } from 'react';
import { Link as NavLink } from 'react-router-dom';
import { SettingsDropdown } from './SettingsDropdown/SettingsDropdown';
import styles from './navigationBar.module.css';
import geneIcon from '/Users/aminbensaad/GeneCruncher-Front/src/assets/logo-white.png';
import profiledefault from '/Users/aminbensaad/GeneCruncher-Front/src/assets/profiledefault.png';


function NavigationBar() {
  const [click, setClick] = useState(false);
  const [dropDown, setDropDown] = useState(false);

  const handleClick = () => setClick(!click);
  const closeMobileMenu = () => setClick(false);

  const onMouseEnter = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(true);
  };
  const onMouseLeave = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(false);
  };

  return (
    <>
      <nav className={styles.navbar}>
        <li className={styles.navItem}>
          <NavLink to="/dashboard">
            <img
              alt="gene-icon"
              src={geneIcon}
              style={{ height: '40px', paddingLeft: '15px' }}
            />
          </NavLink>
        </li>

        <ul
          className={styles.navMenu}
          style={{ left: '50%', margin: 'auto' }}
        >
          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to="/dashboard"
              onClick={closeMobileMenu}
            >
              Dashboard
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to="/documentation"
              onClick={closeMobileMenu}
            >
              Documentation
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to="/help"
              onClick={closeMobileMenu}
            >
              Help
            </NavLink>
          </li>
        </ul>

        <li
          className={styles.navItem}
          onMouseEnter={onMouseEnter}
          onMouseLeave={onMouseLeave}
        >
          <NavLink
            className={styles.profileSettings}
            to="/settings"
            onClick={closeMobileMenu}
          >
            John Doe
            <img
              alt="profiledefault"
              src={profiledefault}
              style={{ height: '40px', paddingLeft: '15px' }}
            />
          </NavLink>
          { dropDown && <SettingsDropdown /> }
        </li>

      </nav>
    </>

  );
}

export default NavigationBar;
