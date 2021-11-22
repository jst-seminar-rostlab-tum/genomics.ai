import React from 'react';

import * as IoIcons from 'react-icons/io';
import * as CGIcons from 'react-icons/cg';

export const SidebarData = [
  {
    id: 0,
    name: 'Project 1',
    path: '/project1',
    // needed in case sub items exist
    iconClosed: <IoIcons.IoIosArrowForward />,
    iconOpened: <IoIcons.IoIosArrowDown />,
    subNav: [
      {
        id: 0,
        name: 'celltype annot. 1',
        path: '/overview/ca1',
        icon: <CGIcons.CgFileDocument />,
      },
      {
        id: 1,
        name: 'celltype annot. 2',
        path: '/overview/ca2',
        icon: <CGIcons.CgFileDocument />,
      },
    ],
  },

  {
    id: 1,
    name: 'Project 2',
    path: '/project2',
    iconClosed: <IoIcons.IoIosArrowForward />,
    iconOpened: <IoIcons.IoIosArrowDown />,
    subNav: [
      {
        name: 'celltype annot. 1',
        path: '/overview/ca1',
        icon: <CGIcons.CgFileDocument />,
      },
      {
        name: 'celltype annot. 2',
        path: '/overview/ca2',
        icon: <CGIcons.CgFileDocument />,
      },
    ],
  },

];

export default SidebarData;
