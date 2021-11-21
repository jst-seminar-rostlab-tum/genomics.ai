import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import styles from './subMenu.module.css';

function SubMenu({ item }) {
  const [subnav, setSubnav] = useState(false);
  const showSubnav = () => setSubnav(!subnav);

  return (
    <>
      <div className={styles.divider} />
      <Link
        className={styles.sidebarLink}
        to={item.path}
        onClick={item.subNav && showSubnav}
      >
        <span className={styles.sidebarLabel}>{ item.name }</span>
        <div>
          { item.subNav && subnav ? item.iconOpened
            : (item.subNav ? item.iconClosed : null) }
        </div>
      </Link>
      {
      subnav && item.subNav.map((item, index) => (
                  <div
                    className={styles.dropdownLink}
                    to={item.path}
                    key={index}
                  >
                    {item.icon}
                    <span className={styles.sidebarLabel}>
                      { item.name }
                    </span>
                  </div>
                ))
            }
    </>
  );
};

export default SubMenu;
