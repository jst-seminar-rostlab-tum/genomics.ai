import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import styles from './subMenu.module.css';

function SubMenu({ item }) {
  const [subnav, setSubnav] = useState(false);
  const showSubnav = () => setSubnav(!subnav);

  function getDropdownStatusIcon() {
    if (item.subNav && subnav) {
      return (
        <>
          {item.iconOpened}
        </>
      );
    }
    return (
      <>
        {item.subNav ? item.iconClosed : null}
      </>
    );
  }

  return (
    <>
      <div className={styles.divider} />
      <Link
        className={styles.sidebarLink}
        to={item.path}
        onClick={item.subNav && showSubnav}
      >
        <span className={styles.sidebarLabel}>{item.name}</span>
        <div>
          {getDropdownStatusIcon()}
        </div>
      </Link>
      {
        subnav && item.subNav.map((subItem) => (
          <Link
            className={styles.dropdownLink}
            to={subItem.path}
            key={subItem.id}
          >
            {subItem.icon}
            <span className={styles.sidebarLabel}>
              {subItem.name}
            </span>
          </Link>
        ))
      }
    </>
  );
}

export default SubMenu;
