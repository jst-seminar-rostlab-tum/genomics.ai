export const SidebarData = [
  {
    id: 0,
    name: 'Project 1',
    path: '/project1',
    // needed in case sub items exist
    subNav: [
      {
        id: 0,
        name: 'celltype annot. 1',
        path: '/overview/ca1',
      },
      {
        id: 1,
        name: 'celltype annot. 2',
        path: '/overview/ca2',
      },
    ],
  },

  {
    id: 1,
    name: 'Project 2',
    path: '/project2',
    subNav: [
      {
        name: 'celltype annot. 1',
        path: '/overview/ca1',
      },
      {
        name: 'celltype annot. 2',
        path: '/overview/ca2',
      },
    ],
  },

];

export default SidebarData;
