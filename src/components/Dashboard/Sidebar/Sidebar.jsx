import React from 'react';
import { SidebarData } from './SidebarData';
import SubMenu from './SubMenu/SubMenu';
import styles from './sidebar.module.css';
import arrowClose from '../../../assets/arrow-close.png';
import geneIcon from '../../../assets/gene.png';

const Sidebar = ({ sidebarShown, toggleSidebar }) => {
  // COLLAPSED
  if (sidebarShown) {
    return (
      <div className={styles.sidebarNavCollapsed}>
        <div className={styles.sidebarWrap}>
          <div>
            <input
              className={styles.projectsBannerCollapsed}
              type="image"
              alt="toggle sidebar"
              src={geneIcon}
              onClick={toggleSidebar}
            />
          </div>
        </div>
      </div>
    );
  }

  // NOT COLLAPSED
  return (
    <div className={styles.sidebarNav}>
      <div className={styles.sidebarWrap}>
        <div style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}>
          <div style={{ paddingBlock: '35px', paddingLeft: '270px' }}>
            <input
              type="image"
              className={styles.toggleButton}
              alt="toggle-icon"
              src={arrowClose}
              onClick={toggleSidebar}
            />
          </div>
        </div>

        <div
          className={styles.projectsBanner}
        >
          <img
            alt="gene-icon"
            src={geneIcon}
            style={{ height: '35px', paddingLeft: '5px', paddingRight: '30px' }}
          />
          Projects
        </div>
        { SidebarData.map((item) => <SubMenu item={item} key={item.id} />)}
      </div>
    </div>
  );
};

export default Sidebar;
